#! /bin/bash
# Run the sc_pileup workflow.
#
# Standalone:   ./run.sh [extra snakemake args]
# As a submodule (installed under bcarc_chromium/modules/sc_pileup), invoke the
# Snakefile with the parent repo's config, mirroring snp_demux's run_mod.sh.

set -e
eval "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "$(hostname)" == *"nci"* ]]; then
    global_profile="--profile $HERE/profiles/global/nci"
    workflow_profile="--workflow-profile $HERE/profiles/workflow"
    module load singularity
    mkdir -p logs/joblogs
else
    global_profile=""
    workflow_profile="--workflow-profile $HERE/profiles/workflow"
fi

snakemake $global_profile $workflow_profile \
    --snakefile "$HERE/workflow/Snakefile" \
    --configfile "${SC_PILEUP_CONFIG:-$HERE/config/config.yaml}" \
    "$@"
