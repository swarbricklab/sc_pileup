# Per-pool informative demux panel.
#
# When config.panel.informative_filter is true, the pileup engines do NOT tally
# at the raw `sites_vcf`. Instead, per capture, this rule builds the sites they
# pile at:
#   1. subset the cohort donor-genotype `sites_vcf` to the pool's donors
#      (captures.csv `donors` column, ';'-separated sample names as they appear
#      in sites_vcf);
#   2. keep only sites that are POLYMORPHIC among those donors -- i.e. sites
#      that can actually distinguish them (`bcftools view -c 1:minor`, which
#      recomputes the allele count on the subset, so it drops sites where the
#      pool is all-hom-ref or all-hom-alt).
#
# Rationale: uninformative sites (all donors share the genotype) contribute ~0
# to the donor likelihood in genotype demux, so dropping them is a near-lossless
# pileup-size reduction. The output file is ALSO intended to be snp_demux's donor
# genotype panel (deps.vcf) -- using the SAME physical file on both sides keeps
# the pileup sites and the demux genotypes identical by construction (they can
# never silently drift and get intersected away).
#
# NB `-c 1:minor` still keeps the rare all-heterozygous site (every donor 0/1),
# which is technically uninformative too; for a common-SNP union panel these are
# a small minority and catching them needs a genotype-level comparison that
# isn't a single bcftools flag, so we accept them.
#
# Leave informative_filter FALSE for somatic / Track-A pileups: there we want
# every site regardless of donor, so engines pile at the raw sites_vcf and this
# rule is never triggered.


rule informative_panel:
    input:
        vcf=config["deps"]["sites_vcf"],
    output:
        vcf=out_dir / "panel/{capture_id}/informative.vcf.gz",
        idx=out_dir / "panel/{capture_id}/informative.vcf.gz.csi",
    params:
        donors=pool_donors,
    log:
        logs / "panel/{capture_id}.log",
    container:
        config["containers"]["demuxafy"]
    shell:
        # ';'-separated in the sheet (CSV-safe) -> ',' for bcftools -s.
        """
        donors=$(echo '{params.donors}' | tr ';' ',')
        bcftools view -s "$donors" -c 1:minor {input.vcf} -Oz -o {output.vcf} 2> {log}
        bcftools index {output.vcf} 2>> {log}
        """
