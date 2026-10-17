#! /bin/bash
# Run the sc_pileup pileup on the bundled PUBLIC test fixture (all three engines),
# locally with singularity, and (with --check) verify it against the DVC-tracked
# golden.
#
# Assumes:
#   * conda env `snakemake_7.32.4`
#   * the fixtures + golden are materialised:  dvc pull   (or dvc checkout)
#   * the two singularity images + genome FASTA exist at the paths in
#     test/config.yaml (NCI-resident; see that file's header)
#
# Usage:
#   ./run_test.sh              # run all three engines into results/test
#   ./run_test.sh --check      # run, then content-compare results/test to test/golden
#   ./run_test.sh -n           # dry run
#   CORES=8 ./run_test.sh      # more local cores
set -euo pipefail

# Make conda usable on a bare (non-login) compute shell, e.g. under `qxub exec`
# in CI: the binary may be on PATH but its shell functions are not set up, so
# `conda activate` errors with "Run 'conda init'". Source the profile if conda is
# absent, then always install the shell hook before activating. Override
# CONDA_PROFILE if your install lives elsewhere.
if ! command -v conda >/dev/null 2>&1; then
    source "${CONDA_PROFILE:-/g/data/a56/conda/miniconda/etc/profile.d/conda.sh}"
fi
eval "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
module load singularity 2>/dev/null || true

# Pull --check out of the args passed through to snakemake.
check=false
sm_args=()
for a in "$@"; do
    if [[ "$a" == "--check" ]]; then check=true; else sm_args+=("$a"); fi
done

snakemake \
    --snakefile "$HERE/workflow/Snakefile" \
    --configfile "$HERE/test/config.yaml" \
    --workflow-profile "$HERE/profiles/workflow" \
    --use-singularity \
    --singularity-args "--bind /g/data/a56,/scratch/a56" \
    --cores "${CORES:-4}" \
    ${sm_args[@]+"${sm_args[@]}"}

if [[ "$check" == true ]]; then
    echo "=== comparing results/test to the golden (test/golden) ==="
    python3 "$HERE/test/compare_golden.py" "$HERE/test/golden" "$HERE/results/test"
fi
