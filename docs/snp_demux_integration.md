# Design — how `snp_demux` picks up where `sc_pileup` leaves off

sc_pileup owns the **pileup** (genotype-free per-cell REF/ALT at a panel);
snp_demux owns **assignment** (donor-aware). This doc defines the interface
contract and the concrete rule changes in snp_demux.

## The interface contract

Per capture `C`, platform `P`, under sc_pileup's results root `R`:

| engine | files | consumed by |
|---|---|---|
| cellsnp | `R/cellsnp/P/C/cellSNP.{base.vcf.gz,samples.tsv,tag.AD.mtx,tag.DP.mtx,tag.OTH.mtx}` | vireo |
| popscle | `R/popscle/P/C/pileup.{cel,var,plp,umi}.gz` | demuxlet **and** freemuxlet |
| souporcell | `R/souporcell/P/C/{ref.mtx,alt.mtx,barcodes.tsv,variants.vcf}` | souporcell clustering |

**Shared invariant (critical):** the panel VCF sc_pileup piles against
(`deps.sites_vcf`) MUST be the *same* donor-keyed genotype panel snp_demux uses.
Then cellSNP `base.vcf.gz` sites line up with the donor genotypes for
`subset_vcf`, and souporcell `variants.vcf` (which carries the donor GT columns,
being a covered subset of that panel) can serve as `--known_genotypes` whose
rows match the `alt.mtx`/`ref.mtx`.

## Per-engine changes in snp_demux

### vireo (cellSNP)

- **Delete** `rule cellsnp_lite`.
- **vireo** input repoints from `rules.cellsnp_lite.output.*` to
  `R/cellsnp/P/C/cellSNP.{base.vcf.gz,samples.tsv,tag.AD.mtx,tag.DP.mtx}`.
- **`subset_vcf` stays** (donor-aware: subsets the panel to the pool's donors at
  the cellSNP sites) — repoint its `base` input to `R/.../cellSNP.base.vcf.gz`.
- `relabel_vireo`, `reheader_vireo_vcf`, `summarise_ambient` — unchanged.

### demuxlet + freemuxlet (popscle)

- **Delete** `split_bam_by_chromosome`, `split_vcf_by_chromosome`, both
  `dsc_pileup_*_by_chromosome`, and every `merge_*` / `wrangle_*` rule (and the
  `merge_cel/merge_var/wrangle_plp/wrangle_umi` scripts).
- **demuxlet** input: pileup prefix `R/popscle/P/C/pileup` (`.plp/.cel/.var/.umi.gz`)
  + `known_genotypes` VCF (donor subset — stays). The old `--sm-list` on the
  pileup is gone; it was inert to pileup content, and demuxlet still restricts
  via `--vcf` at assignment.
- **freemuxlet** input: the *same* shared pileup (this is where the de-dup pays
  off — one pileup, both tools).
- `demuxlet_summary`, `freemuxlet_genotyping`, summaries — unchanged.

### souporcell

Replace the monolithic `rule souporcell` (remap+vartrix+cluster+troublet+
consensus) with three lean rules that resume from sc_pileup's matrices:

```
# 1. cluster (donor-aware: k = #donors; known genotypes guide clustering)
souporcell -k {n_donors} -a alt.mtx -r ref.mtx -b barcodes.tsv \
    --restarts {N} --min_ref {min_ref} --min_alt {min_alt} \
    --known_genotypes variants.vcf --known_genotypes_sample_names {pool_donors} \
    > clusters_tmp.tsv
# 2. doublets
troublet --alts alt.mtx --refs ref.mtx --clusters clusters_tmp.tsv > clusters.tsv
# 3. consensus genotypes + ambient
consensus.py -c clusters.tsv -a alt.mtx -r ref.mtx --soup_out ambient_rna.txt \
    -v variants.vcf --output_dir . --vcf_out cluster_genotypes.vcf
```

- inputs: `R/souporcell/P/C/{alt.mtx,ref.mtx,barcodes.tsv,variants.vcf}`.
- `n_donors` = `get_number_of_samples_in_capture_pool`; `pool_donors` =
  `get_donors_in_capture` (donor-aware helpers — stay).
- **`--known_genotypes` must be sc_pileup's `variants.vcf`** (the covered subset)
  so its rows match the matrix rows; restrict to the pool via
  `--known_genotypes_sample_names`. This *replaces* the old `known_genotypes`
  subset rule for souporcell (that VCF had all panel sites and would not align
  with the covered matrices).
- Outputs (`clusters.tsv`, `cluster_genotypes.vcf`, `ambient_rna.txt`, plus
  `alt.mtx`/`ref.mtx` now sourced from sc_pileup) match what the current stage
  keeps, so `souporcell_genotyping.smk` (vireosnp cluster→donor) is unchanged.

### demuxafy consensus, aggregate, report

Unchanged — they consume the per-method assignment outputs, which keep the same
paths and formats.

## What stays in snp_demux (all donor-aware)

`subset_vcf`, `known_genotypes`/`extract_genotyped_donors` (for demuxlet's donor
list), `get_number_of_samples_in_capture_pool` (k), `get_donors_in_capture`, and
everything from assignment → relabel → consensus → aggregate → report.

## Wiring in bcarc_chromium

- New **`sc_pileup` DVC stage** producing `data/snp/pileup/{cellsnp,popscle,souporcell}/P/C/...`
  (tracked outs). Runs `modules/sc_pileup/run.sh`.
- **`snp_demux` stage** gains a dep on `data/snp/pileup`; its rules read from a
  configured `deps.pileup_root: data/snp/pileup`. Stage order: `sc_pileup` →
  `snp_demux`.
- **config**: sc_pileup `deps.sites_vcf` = `data/snp/vcf/donor_keyed.vcf.gz`
  (the SAME panel snp_demux uses). snp_demux `captures.csv` needs a `platform`
  column (default `short`) to locate the pileups.
- Both stages are currently `frozen: true` — coordinate an un-freeze/rebuild.

## Migration & validation

- Add a `pileup_source: internal | sc_pileup` switch in snp_demux so it can run
  either way during transition, enabling an **A/B check**: demux assignments
  should be **identical** (same pileups feed the same assigners).
  - cellSNP matrices are byte-for-byte reproducible (same params/BAM/VCF) →
    vireo identical.
  - popscle shared pileup → demuxlet/freemuxlet identical (`--sm-list` inert).
  - souporcell prep is faithful → clustering equivalent, **modulo** souporcell's
    intrinsic `--restarts` stochasticity (cluster *labels* may permute; the
    vireosnp cluster→donor step resolves that, so final donor assignments match).
    Note this nondeterminism already exists in the monolith; the split adds none.
- Success metric: per-cell assignment concordance ≈ 100% vs the pre-split stage,
  especially among epithelial cells.

## Gotchas

1. **Panel VCF identity** between sc_pileup and snp_demux is mandatory.
2. **souporcell `--known_genotypes` = `variants.vcf`** (covered subset), not the
   full-panel donor subset — else sites won't match the matrix rows.
3. **Platform dimension**: snp_demux must know each capture's platform to find
   the right pileup dir (default `short` until long-read lands).
4. **Frozen DVC stages**: rebuild both in one coordinated pass; `dvc repro`
   deletes then rebuilds `data/snp/demux`.
5. **barcodes consistency**: every engine's barcodes derive from the same
   filtered `barcodes.tsv.gz`, so cell axes stay aligned across cellSNP
   (`samples.tsv`), souporcell (`barcodes.tsv`), and popscle (`cel`).
