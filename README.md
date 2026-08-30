# sc_pileup

A small, reusable Snakemake workflow that turns single-cell BAMs into **per-cell
REF/ALT count matrices** at a configurable set of sites. It factors the
expensive, deterministic *pileup* out of the *assignment* step so that:

- **genotype demux** (vireo / demuxlet / freemuxlet in `snp_demux`) can be
  re-run and re-tuned cheaply, resuming from saved pileups instead of re-reading
  BAMs;
- the same matrices feed **malignant-cell / clonal** analyses — somatic-variant
  detection (Track A), allelic imbalance in CNV/LOH regions (Track B), and
  `cardelino` — when pointed at a somatic / CNV-region site VCF.

Design rationale and the wider plan (WGS+microarray panel, demux improvements,
cardelino) live in [`docs/design_pileup_stage.md`](docs/design_pileup_stage.md).

## The abstraction

```
pileup(bam, barcodes, platform, sites_vcf) -> per-cell REF/ALT counts at sites
```

The pileup is **genotype-free**: it tallies every barcode at every site,
irrespective of donor. Donor subsetting / assignment happens downstream. So the
only sample sheet is `config/captures.csv` (one row per capture); no per-donor
genotype wrangling lives here.

## Engines

| engine | output | consumers | platforms |
|---|---|---|---|
| `cellsnp` | cellSNP-lite MEX: `cellSNP.base.vcf.gz`, `cellSNP.samples.tsv`, `cellSNP.tag.{AD,DP,OTH}.mtx` | vireo, Track B, cardelino | short, long* |
| `popscle` | `dsc-pileup`: `pileup.{cel,var,plp,umi}.gz` | demuxlet **and** freemuxlet (shared) | short |

Outputs are written under `results/<engine>/<platform>/<capture_id>/` so
short- and long-read pileups of the same capture coexist and can be merged
per cell downstream.

\* **Long-read (PacBio) note.** Full-length HiFi reads cover the whole transcript
body, not just the 3′ end, so many more exonic somatic sites become observable —
this is the main reason the malignant-detection direction is worth pursuing. The
`platform` dimension is wired through the config and output layout, but the
long-read pileup parameters (splice-aware handling, UMI semantics, SNV-preferring
site QC) are **placeholders to be tuned against real PacBio single-cell BAMs**
(`params.cellsnp_lite.platform_args.long`). `popscle` is short-read only for now.

## Provenance

The pileup rules are heavily plagiarised from
[`swarbricklab/snp_demux`](https://github.com/swarbricklab/snp_demux)
(`rules/vireo.smk` cellSNP, `rules/demuxlet.smk`/`freemuxlet.smk` dsc-pileup,
`rules/split.smk`, and `scripts/{merge_cel,merge_var,wrangle_plp,wrangle_umi}.py`),
then reviewed and generalised. Key changes: outputs are no longer `temp()`; the
site VCF is generic; the redundant demuxlet/freemuxlet pileups are merged into
one; genotype-free (no `--sm-list`); platform-aware paths. See the header
comments in each rule file.

## Quickstart

1. Edit `config/config.yaml`: point `deps.sites_vcf` at a bgzipped+indexed VCF
   (GRCh38, matching the BAMs), choose `engines`, and set the cellsnp `minMAF`/
   `minCOUNT` for your purpose (germline panel: `0.01`/`10`; targeted somatic:
   `0`/`1`).
2. Fill `config/captures.csv`: `capture_id,platform,bam,barcodes`.
3. Provide the container at `resources/Demuxafy.sif` (same image `snp_demux`
   uses).
4. Dry run, then run:
   ```bash
   ./run.sh -n          # dry run
   ./run.sh -c 32       # local
   ```

## Status

Experimental. Intended to be installed as a git submodule under
`bcarc_chromium/modules/sc_pileup` on an experimental branch once tested.
