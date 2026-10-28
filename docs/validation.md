# Validation

sc_pileup has two kinds of validation:

1. **CI fixture + golden (determinism).** `run_test.sh --check` runs all three
   engines on the DVC-tracked 10%-downsampled OneK1K fixture and content-compares to
   the frozen golden. This proves the pileup is byte-reproducible. See
   `.tests/integration/make_test_fixtures.sh` for fixture provenance.

2. **End-to-end correctness (consistency with the Demuxafy reference).** sc_pileup's
   per-cell pileups are the input to `snp_demux`'s demultiplexing. On the **full**
   OneK1K dataset, running `sc_pileup → snp_demux` reproduces the per-method demux
   numbers that Neavin et al. publish for Demuxafy (vireo/freemuxlet/souporcell).
   Because the demux calls depend directly on the pileup's REF/ALT counts, that
   agreement validates the pileup itself end-to-end. The full write-up (datasets,
   per-method comparison tables, the matched-cellSNP apples-to-apples for vireo) lives
   in **snp_demux: `docs/validation/demuxafy_consistency.md`**.

## Full-depth run configs (for the correctness validation)

The heavy, one-off full-depth pileup is not part of CI (the 39 GB input is external,
not DVC-tracked). Its configs are committed for reproducibility:

- `.tests/integration/config.full.yaml` — all three engines on the full BAM; cellSNP kept lenient
  (`minMAF 0.0 / minCOUNT 1`) to retain every covered site.
- `.tests/integration/config.full_neavin.yaml` — cellSNP only, at **Neavin's thresholds**
  (`minMAF 0.1 / minCOUNT 20`), used for the apples-to-apples vireo comparison.

Both point at a full fixture built by the same transform as the CI fixture
(`.tests/integration/make_test_fixtures.sh`) but without the 10% read subsample: chr-renamed BAM +
the 14-donor genotyped sites VCF. Run on PBS with `--use-singularity` (containers
pre-cached in the shared prefix; compute nodes have no internet).
