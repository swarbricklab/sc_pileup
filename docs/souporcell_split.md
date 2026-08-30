# Design — splitting souporcell between `sc_pileup` and `snp_demux`

## Motivation

souporcell's `souporcell_pipeline.py` is a monolith of six stages:

1. **remap** — extract reads to fastq (`renamer.py`), realign with minimap2
   (splice, ref-bias-reducing), restore CB/UB tags (`retag.py`), sort+index
2. **variant calling** — freebayes, *or* skipped when a variant VCF is supplied
   (`--common_variants` / `--known_genotypes`), reduced to the **covered** subset
   by a samtools-depth × bedtools-intersect pass
3. **allele counting** — `vartrix` → `ref.mtx` / `alt.mtx` (cell × variant)
4. **clustering** — the `souporcell` Rust binary (`-k` = #donors)
5. **doublets** — `troublet`
6. **consensus** — `consensus.py` (genotypes + ambient RNA)

Stages 1–3 are expensive, deterministic, BAM-reading, and **donor-count-agnostic**
— exactly a pileup's job. Stages 4–6 are fast, depend on `k`/known genotypes, and
are **donor-aware** — demux's job. Today snp_demux runs all six as one frozen,
~4 h/capture rule, so a demux re-tune re-does the remap. This mirrors the
cellSNP/popscle split we already made: **pileup once, resume demux cheaply.**

## The cut point: `ref.mtx` / `alt.mtx`

souporcell documents running the clustering binary directly on precomputed
matrices:

```
souporcell -a alt.mtx -r ref.mtx -b barcodes.tsv -k <num_clusters> -t 8 > clusters.tsv
```

So the boundary is clean:

| stage | module | notes |
|---|---|---|
| renamer → minimap2 → retag → sort | **sc_pileup** | the remap (ref-bias correction) |
| covered-variant filter | **sc_pileup** | samtools depth × bedtools intersect on the *remapped* bam |
| vartrix → `ref.mtx`/`alt.mtx` | **sc_pileup** | + `barcodes.tsv`, `variants.vcf` kept alongside |
| `souporcell` clustering (`-k`) | **snp_demux** | needs #donors |
| troublet, consensus | **snp_demux** | doublets, cluster→donor, ambient |

sc_pileup runs the prep **genotype-free**: vartrix at the shared `sites_vcf`
panel (via `--common_variants` semantics), *not* `--known_genotypes`. The number
of clusters and the optional cluster→donor mapping (`--known_genotypes`) stay in
snp_demux's clustering step. This keeps the pileup donor-agnostic, consistent
with cellSNP/popscle.

### Bonus: this fixes snp_demux issue #35 for souporcell

snp_demux keeps `alt.mtx`/`ref.mtx` but `temp()`s their `barcodes.tsv` and the
variant VCF that index them, leaving orphaned matrices. Here sc_pileup emits all
four together (`ref.mtx`, `alt.mtx`, `barcodes.tsv`, `variants.vcf`) as one
coherent, non-`temp()` product.

## Faithful commands (from souporcell_pipeline.py, v2.0-66 ≡ 2.5 for these stages)

- **renamer**: `renamer.py --bam BAM --barcodes BC --out reads.fq --no_umi False --umi_tag UB --cell_tag CB`
- **minimap2**: `minimap2 -ax splice -t N -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no FASTA reads.fq > minimap.sam`
- **retag**: `retag.py --sam minimap.sam --no_umi False --umi_tag UB --cell_tag CB --out tagged.bam` → `samtools sort`+`index`
- **covered variants** (on the *remapped* bam): `samtools depth` → keep `depth ≥ min_cov` → `bedtools merge` → `bedtools intersect -wa -a sites_vcf -b covered.bed` (prepend the VCF header)
- **vartrix**: `vartrix --mapq 30 -b remapped_sorted.bam -c barcodes.tsv --scoring-method coverage --threads N --ref-matrix ref.mtx --out-matrix alt.mtx -v covered.vcf --fasta FASTA --umi`

## Chromosome-split optimisation — partial

- **minimap2 remap**: **not** cleanly chrom-splittable — remapping's purpose is to
  let mis-mapped reads relocate genome-wide; splitting reads by original
  chromosome would miss cross-chromosome relocations. minimap2 is already
  internally multithreaded.
- **covered-variant depth + vartrix**: per-locus, so chrom-splittable (souporcell
  itself recommends *"1 job per chromosome and merge after"* for variant work).
  Merging the per-chrom `ref.mtx`/`alt.mtx` needs MatrixMarket row-reindexing —
  the same trick the popscle `merge_*` scripts do for SNP_ID. **Deferred to v1.1.**

v1 runs whole-BAM. The chrom-split lever helps only the post-remap half.

## Version: pin 2.x now, evaluate souporcell3 separately

Upstream has bumped to **souporcell3** (Bioinformatics 2024). Do **not** bundle
that upgrade with this refactor:

- the split must be **behaviour-preserving**, so it pins the same toolchain
  snp_demux clusters with — souporcell/vartrix/minimap2 from `Demuxafy.sif`
  (vartrix 1.1.3, minimap2 2.7), whose prep commands match `cumulusprod/souporcell:2.5`;
- souporcell3 may change matrix/VCF formats and cluster output, which would
  ripple into snp_demux's clustering, `vireosnp` annotation, Demuxafy consensus,
  and `geno_concordance.py` — and bcarc_chromium's demux results already exist on
  2.5. Evaluate v3 as its own experiment, gated on a demux-concordance check
  (does v3 reproduce v2.5 assignments?).

## How snp_demux would consume this (future, in that repo)

Replace the monolithic souporcell rule with a clustering rule that resumes from
sc_pileup's outputs:

```
souporcell -a alt.mtx -r ref.mtx -b barcodes.tsv -k <n_donors> -t N > clusters_tmp.tsv
troublet --alts alt.mtx --refs ref.mtx --clusters clusters_tmp.tsv > clusters.tsv
consensus.py -c clusters.tsv -a alt.mtx -r ref.mtx --soup_out ambient.txt \
    -v variants.vcf --output_dir . --vcf_out cluster_genotypes.vcf
```

Everything downstream (`cluster_genotypes.vcf`, vireosnp cluster→donor, Demuxafy
consensus) is unchanged. Only the boundary moves.

## Optional future optimisation: reuse cellSNP instead of a separate vartrix

With `--skip_remap`, souporcell-prep collapses to "vartrix at panel sites on the
original BAM" — conceptually what cellSNP already produces (alt = AD, ref =
DP−AD). One could then feed clustering from cellSNP matrices (format-convert to
MatrixMarket) and run *no* separate souporcell pileup. Trade-off: loses the
remap-based ref-bias correction. Benchmark before adopting; v1 keeps remap.
