# Somatic site-VCF checks (panel_kind: somatic only).
#
# cellSNP-lite piles up biallelic SNVs, so for a somatic (Track A) site VCF we
# count the records that are not (indels, MNVs, multiallelic, symbolic) before
# piling up. `somatic.non_snv: warn` (default) reports them in the log and in
# sites/<capture>/site_types.tsv; `drop` also writes a biallelic-SNV-only subset
# that the engines then pile at (see sites_vcf_input). Germline captures never
# reach these rules, so their DAG is unchanged.


rule somatic_site_types:
    input:
        vcf=raw_sites_vcf,
        script=str(Path(workflow.basedir) / "scripts/check_site_types.sh"),
    output:
        report=out_dir / "sites/{capture_id}/site_types.tsv",
    log:
        logs / "sites/{capture_id}/site_types.log",
    container:
        config["containers"]["bcftools"]
    params:
        mode=_somatic_non_snv,
    shell:
        """
        bash {input.script} {input.vcf} {params.mode} {wildcards.capture_id} \
            {output.report} >{log} 2>&1
        # surface the counts (and any warning) in the run's console output too
        cat {log} >&2
        """


rule somatic_snv_sites:
    input:
        vcf=raw_sites_vcf,
        report=rules.somatic_site_types.output.report,
        script=str(Path(workflow.basedir) / "scripts/check_site_types.sh"),
    output:
        vcf=out_dir / "sites/{capture_id}/biallelic_snv.vcf.gz",
        idx=out_dir / "sites/{capture_id}/biallelic_snv.vcf.gz.csi",
    log:
        logs / "sites/{capture_id}/biallelic_snv.log",
    container:
        config["containers"]["bcftools"]
    shell:
        """
        d=$(dirname {output.vcf})
        bash {input.script} {input.vcf} drop {wildcards.capture_id} \
            $d/.drop_report.tsv {output.vcf} >{log} 2>&1
        rm -f $d/.drop_report.tsv
        """
