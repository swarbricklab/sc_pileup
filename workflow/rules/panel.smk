# Optional in-workflow per-pool panel build (STANDALONE convenience).
#
# This rule is a thin wrapper around tools/build_informative_panel.sh so that a
# user running sc_pileup on its own (just a cohort/project donor-genotype VCF, no
# orchestrator) can still get per-pool demux panels. It is triggered only when
# config.panel.build_in_workflow is true; engines then pile at the built panel
# (see sites_vcf_input, mode 1).
#
# In production the SAME tool is instead called by an upstream orchestration DVC
# stage (e.g. bcarc_chromium build_informative_panels) in batch across all
# captures, and sc_pileup consumes those panels via a per-capture `sites_vcf`
# column (mode 2) -- so there is ONE implementation, two entry points.
#
# Steps (all optional, config-driven; folded into one bcftools view by the tool):
#   panel.common_sites -> restrict to a common-SNP backbone (harmonise sources)
#   donors             -> subset to the pool (captures.csv `donors` column)
#   panel.informative  -> keep only sites polymorphic among those donors
# The output doubles as snp_demux's donor genotype panel (same file => pileup
# sites and demux genotypes identical by construction).


rule informative_panel:
    input:
        vcf=config["deps"]["sites_vcf"],
        tool=str(tools_dir / "build_informative_panel.sh"),
    output:
        vcf=out_dir / "panel/{capture_id}/informative.vcf.gz",
        idx=out_dir / "panel/{capture_id}/informative.vcf.gz.csi",
    params:
        donors=pool_donors,
        out_root=str(out_dir / "panel"),
        common_arg=panel_common_arg,
        informative_arg=panel_informative_arg,
    log:
        logs / "panel/{capture_id}.log",
    container:
        config["containers"]["bcftools"]
    shell:
        """
        bash {input.tool} \
            --vcf {input.vcf} \
            --out-dir {params.out_root} \
            --capture {wildcards.capture_id} \
            --donors '{params.donors}' \
            {params.common_arg} {params.informative_arg} \
            > {log} 2>&1
        """
