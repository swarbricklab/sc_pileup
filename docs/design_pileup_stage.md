# Design 1 — A standalone `pileup` stage feeding `snp_demux`

## Goal

Decouple the expensive, deterministic **pileup** (reading the BAM to get per-cell
REF/ALT counts at a set of sites) from the cheap, tunable **assignment** step
(vireo/souporcell/demuxlet/freemuxlet). Benefits:

1. **Fast demux reruns** — re-tuning thresholds, panels, or consensus logic no
   longer re-reads the BAMs (currently ~4 h/capture, souporcell the long pole).
2. **Reuse** — the same per-cell count matrices feed the malignant-cell work
   (Track A somatic detection, Track B allelic imbalance, cardelino).
3. **Dedup** — `demuxlet` and `freemuxlet` currently each run their own popscle
   pileup over the same BAM; a shared pileup halves that.

This complements snp_demux issue #35 (stop `temp()`-ing the pileups): that issue
keeps the *intermediates*; this doc promotes them to a *first-class stage*.

## The abstraction

```
pileup(BAM, barcodes, site_VCF) -> per-cell REF/ALT counts at site_VCF sites
```

Genotype-free: cellSNP-lite / popscle `dsc-pileup` tally *every* barcode at
*every* panel site irrespective of donor. So **one pileup per (capture, panel)
serves all donors and all downstream demux/malignancy tools** — donor subsetting
happens later, on the counts.

## Reality check: three consumers, two pileup formats

| Consumer | Pileup it needs | Cleanly resumable? |
|---|---|---|
| vireo | cellSNP-lite MEX (`base.vcf.gz` + `AD/DP/OTH.mtx`) | **Yes** — vireo already reads a cellSNP dir |
| demuxlet | popscle `dsc-pileup` (`.plp/.cel/.var/.umi`) | **Yes** — popscle is natively 2-step |
| freemuxlet | popscle `dsc-pileup` (same files) | **Yes** — can share demuxlet's pileup |
| souporcell | internal minimap2 remap → vartrix | **No** — monolithic; cache `alt.mtx`/`ref.mtx` |

So the `pileup` stage produces **two artifacts per capture**:

- `pileup/cellsnp/{capture}/` — cellSNP-lite MEX → vireo (+ reusable for Track B / cardelino)
- `pileup/popscle/{capture}/` — one `dsc-pileup` → demuxlet **and** freemuxlet

souporcell stays as-is; treat its `alt.mtx`/`ref.mtx` (+ the barcode/variant
indices from issue #35) as its cached "pileup".

## Proposed DAG

```
                 ┌─────────────────────────────────────────────┐
  panel VCF ───▶ │  pileup stage (DVC-tracked, saved)           │
  BAM        ───▶│   • cellsnp_lite  → pileup/cellsnp/{cap}/     │
  barcodes   ───▶│   • dsc_pileup    → pileup/popscle/{cap}/     │
                 └─────────────────────────────────────────────┘
                        │                 │
                        ▼                 ▼
                 snp_demux: vireo   snp_demux: demuxlet + freemuxlet
                 (resume)           (resume, shared pileup)
                        │                 │
                        └──── demuxafy consensus ────┘
   (souporcell runs its own remap+vartrix branch in parallel, matrices cached)
```

In DVC terms: a new `pileup` stage upstream of `snp_demux`; `snp_demux`'s `deps`
gain `data/snp/pileup/**` and its rules point at those instead of recomputing.
Add `data/snp/pileup` to the consuming repo's `dvc.yaml` `outs` (note from #35:
`data/snp/demux/cellSNP` is *not* currently tracked — the new location must be).

## Rule sketch (Snakemake)

```python
# rules/pileup.smk  (new module stage, or a pre-stage in snp_demux)
rule cellsnp_pileup:
    input:
        bam=lookup(captures, "bam"),
        barcodes=lookup(captures, "barcodes"),
        vcf=config['deps']['panel_vcf'],          # see "which VCF" below
    output:
        base=pileup_dir/"cellsnp/{capture_id}/cellSNP.base.vcf.gz",
        samples=pileup_dir/"cellsnp/{capture_id}/cellSNP.samples.tsv",
        ad=pileup_dir/"cellsnp/{capture_id}/cellSNP.tag.AD.mtx",
        dp=pileup_dir/"cellsnp/{capture_id}/cellSNP.tag.DP.mtx",
        oth=pileup_dir/"cellsnp/{capture_id}/cellSNP.tag.OTH.mtx",
    params:
        minMAF=config['params']['cellsnp_lite']['minMAF'],
        minCOUNT=config['params']['cellsnp_lite']['minCOUNT'],
    container: config["containers"]["demuxafy"]
    shell:
        "cellsnp-lite --samFile {input.bam} --regionsVCF {input.vcf} "
        "--barcodeFile {input.barcodes} --outDir $(dirname {output.base}) "
        "--nproc {resources.threads} --minMAF {params.minMAF} "
        "--minCOUNT {params.minCOUNT} --gzip"

rule dsc_pileup:                                   # feeds BOTH demuxlet & freemuxlet
    input:
        bam=lookup(captures, "bam"),
        barcodes=lookup(captures, "barcodes"),
        vcf=config['deps']['panel_vcf'],
    output:
        multiext(str(pileup_dir/"popscle/{capture_id}/pileup"),
                 ".plp.gz", ".cel.gz", ".var.gz", ".umi.gz"),
    container: config["containers"]["demuxafy"]
    shell:
        "popscle dsc-pileup --sam {input.bam} --vcf {input.vcf} "
        "--group-list {input.barcodes} --out $(dirname {output[0]})/pileup"
```

vireo then consumes `pileup/cellsnp/{capture}/` directly; demuxlet/freemuxlet
take `--plp pileup/popscle/{capture}/pileup`. (popscle's per-chromosome split
+ merge in the current `split.smk` collapses into one `dsc-pileup` call.)

## Which VCF to pile against

The pileup stage should accept the panel VCF as a **config parameter**, and we
run it with different panels for different purposes:

| Purpose | Panel VCF | Notes |
|---|---|---|
| Donor demux (default) | germline common-SNP panel (`donor_keyed.vcf.gz`, upgraded — see below) | dense, shared across donors |
| Track A (somatic detection) | per-donor somatic VCF (`*.purple.somatic.vcf.gz`) | `--minMAF 0 --minCOUNT 1`; per-donor, so per-capture-donor pileups |
| Track B (allelic imbalance) | the germline demux pileup, later **subset** to het sites in PURPLE LOH/gain segments | free by-product of demux pileup |

Track B needs no new pileup — it reuses the germline cellSNP matrices and
filters sites downstream. Track A does need its own (somatic sites aren't in the
germline panel). Same rule, different `panel_vcf`.

## Upgrading the demux panel: WGS + microarray

Today the demux panel is the **imputed microarray** VCF, reheadered to bare
donor IDs:

- `data/snp/imputation/imputed_filtered.hg38.vcf.gz` — multi-sample, imputed,
  dense common SNPs, sample names `SNP_<donor>_<tissue>_S<n>`; covers all
  arrayed donors. Reheadered → `data/snp/vcf/donor_keyed.vcf.gz`.

WGS germline gives something better *for the donors that have it*:

- `bcarc_wgs/data/tumor-normal/<run>/<NORMALLIB>_dragen_germline/<NORMALLIB>.hard-filtered.vcf.gz`
  — single-sample, exact (not imputed), includes rare variants, sample name =
  library ID `L#######`. Covers only WGS'd donors (31 of the neochemo_sc set).

### Why combine

- WGS is **exact and denser** (real genotypes incl. rare variants) → more
  informative sites → cleaner donor separation, especially where imputation is
  weak.
- Microarray is **universal** — the only genotype source for the ~25 neochemo_sc
  donors without WGS. Dropping them isn't an option in a pooled capture.

### Combining strategy (fixed-panel, best-source fill)

Keep a **single fixed site panel** (demux tools behave best when all donors are
genotyped on the same sites), and fill each donor's genotypes from the best
available source:

1. **Define the panel** = the microarray/imputed site list (already genome-wide,
   ~millions of common SNPs, consistent across donors). Optionally intersect
   with 1000G common SNPs for robustness.
2. **Per donor, choose the source:**
   - WGS available → take WGS germline genotypes at the panel sites
     (`bcftools view -R panel.bed <NORMALLIB>.hard-filtered.vcf.gz`), keeping
     `0/1`,`1/1`,`0/0`.
   - else → microarray/imputed genotypes at the panel sites.
3. **Harmonise identity & alleles:**
   - Crosswalk `library → donor_id` via `bcarc_wgs/neochemo.csv` (+ `index.csv`);
     microarray `gussid → donor_id` via `config/snp_demux/donor_to_gussid.tsv`.
   - Reheader every source to bare `donor_id`.
   - Align REF/ALT and normalise (`bcftools norm -f genome.fa`) so WGS and array
     rows are mergeable; drop sites where they can't be reconciled.
4. **Merge** into one multi-sample `donor_keyed.vcf.gz` keyed by donor_id
   (`bcftools merge`); where a donor is missing a site, leave `./.` — vireo and
   souporcell tolerate missing genotypes; demuxlet needs the site present for at
   least the pool's donors.

### Practical notes / gotchas

- **Don't upgrade only some donors' sites** — if WGS donors get extra sites that
  array donors lack, those sites become uninformative for the pool (missing in
  half the donors). Sticking to the shared panel and only swapping the *genotype
  source* avoids this. (A "WGS-only extra sites" panel can be a *separate*,
  optional high-resolution panel used when a pool is all-WGS.)
- **Sanity-check with the existing concordance work** — this repo's
  `geno_concordance.py` already compares imputed vs inferred genotypes; extend it
  to WGS-vs-array per donor to catch swaps/mislabels *before* they poison demux
  (the README already documents a genotype-swap chain).
- **Reheader once, upstream** — building the merged panel belongs in the panel-
  prep step (analogous to the current `snp_reheader_vcf` stage), not in the
  pileup rule. The pileup stage just consumes `config['deps']['panel_vcf']`.

## Migration / rollout

1. Land issue #35 (stop deleting pileups) — minimal change.
2. Add `rules/pileup.smk` + a `pileup` DVC stage; repoint vireo and
   demuxlet/freemuxlet at the saved pileups; collapse the redundant
   demuxlet/freemuxlet pileup rules.
3. Build the merged WGS+array panel in a panel-prep step; switch
   `deps.panel_vcf` to it; re-demux and compare assignment rates (esp. among
   epithelial cells) against the current microarray-only panel.
4. Add a second pileup invocation with the per-donor somatic VCF for Track A /
   cardelino (Design 3).
```
