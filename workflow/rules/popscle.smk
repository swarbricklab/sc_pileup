# popscle dsc-pileup — feeds BOTH demuxlet and freemuxlet downstream.
#
# Plagiarised from snp_demux (demuxlet.smk / freemuxlet.smk / split.smk) and
# de-duplicated: snp_demux ran a separate dsc-pileup for demuxlet and again for
# freemuxlet; here ONE pileup per capture serves both.
#
# Reviewed changes vs snp_demux:
#   * genotype-free: no `--sm-list` (that was a demuxlet-only optimisation);
#     the pileup tallies all barcodes at all sites, donor subsetting is
#     downstream;
#   * site VCF is generic `sites_vcf`;
#   * final merged pileup.{cel,var,plp,umi}.gz are NOT temp() (the product);
#     per-chromosome chunks and index maps stay temp();
#   * paths carry {platform} (popscle currently short-read only — see config
#     params.popscle.platforms).
#
# NOTE: chr1..chr22 only (no chrX/Y), matching snp_demux. Revisit if sex
# chromosomes are needed for a given panel.

# NB: plain concatenation, NOT f-strings — Snakemake's source reconstruction
# inserts spaces around `{...}` in f-strings inside .smk files (it would turn
# f"chr{i}" into " chr 5 ").
chromosomes = ["chr" + str(i) for i in range(1, 23)]


rule split_bam_by_chromosome:
    input:
        bam=lookup("bam"),
    output:
        bam=temp(out_dir / "popscle/{platform}/{capture_id}/split/{chrom}.bam"),
        index=temp(out_dir / "popscle/{platform}/{capture_id}/split/{chrom}.bam.bai"),
    log:
        logs / "popscle/{platform}/{capture_id}/split_bam_{chrom}.log",
    container:
        config["containers"]["samtools"]
    shell:
        """
        samtools view -b {input.bam} {wildcards.chrom} >{output.bam} 2>{log}
        samtools index {output.bam} 2>>{log}
        """


rule split_vcf_by_chromosome:
    input:
        vcf=sites_vcf_input,
    output:
        vcf=temp(out_dir / "popscle/{platform}/{capture_id}/split/{chrom}.vcf.gz"),
        index=temp(out_dir / "popscle/{platform}/{capture_id}/split/{chrom}.vcf.gz.csi"),
    log:
        logs / "popscle/{platform}/{capture_id}/split_vcf_{chrom}.log",
    container:
        config["containers"]["bcftools"]
    shell:
        """
        bcftools view -r {wildcards.chrom} {input.vcf} -Oz -o {output.vcf} 2>{log}
        bcftools index {output.vcf} 2>>{log}
        """


rule dsc_pileup_by_chromosome:
    input:
        bam=rules.split_bam_by_chromosome.output.bam,
        vcf=rules.split_vcf_by_chromosome.output.vcf,
        barcodes=lookup("barcodes"),
    output:
        cel=temp(out_dir / "popscle/{platform}/{capture_id}/by_chrom/{chrom}.cel.gz"),
        plp=temp(out_dir / "popscle/{platform}/{capture_id}/by_chrom/{chrom}.plp.gz"),
        umi=temp(out_dir / "popscle/{platform}/{capture_id}/by_chrom/{chrom}.umi.gz"),
        var=temp(out_dir / "popscle/{platform}/{capture_id}/by_chrom/{chrom}.var.gz"),
    log:
        logs / "popscle/{platform}/{capture_id}/dsc_pileup_{chrom}.log",
    container:
        config["containers"]["popscle"]
    params:
        min_mq=config["params"]["popscle"]["min_mq"],
        min_bq=config["params"]["popscle"]["min_bq"],
        cap_bq=config["params"]["popscle"]["cap_bq"],
        excl_flag=config["params"]["popscle"]["excl_flag"],
        extra=config["params"]["popscle"]["extra"],
    shell:
        """
        cel={output.cel}
        out_prefix=${{cel%.cel.gz}}
        popscle dsc-pileup \
            --sam {input.bam} \
            --vcf {input.vcf} \
            --group-list {input.barcodes} \
            --min-MQ {params.min_mq} \
            --min-BQ {params.min_bq} \
            --cap-BQ {params.cap_bq} \
            --excl-flag {params.excl_flag} {params.extra} \
            --out $out_prefix \
            >{log} 2>&1
        """


rule merge_cel:
    input:
        cel=expand(
            out_dir / "popscle/{{platform}}/{{capture_id}}/by_chrom/{chrom}.cel.gz",
            chrom=chromosomes,
        ),
    output:
        cel=out_dir / "popscle/{platform}/{capture_id}/pileup.cel.gz",
        sample_map=temp(
            out_dir / "popscle/{platform}/{capture_id}/full_sample_map.tsv.gz"
        ),
    log:
        logs / "popscle/{platform}/{capture_id}/merge_cel.log",
    script:
        "../scripts/merge_cel.py"


rule merge_var:
    input:
        var=expand(
            out_dir / "popscle/{{platform}}/{{capture_id}}/by_chrom/{chrom}.var.gz",
            chrom=chromosomes,
        ),
    output:
        var=out_dir / "popscle/{platform}/{capture_id}/pileup.var.gz",
        snp_map=temp(out_dir / "popscle/{platform}/{capture_id}/snp_map.tsv.gz"),
    log:
        logs / "popscle/{platform}/{capture_id}/merge_var.log",
    script:
        "../scripts/merge_var.py"


rule wrangle_plp_by_chromosome:
    input:
        plp=out_dir / "popscle/{platform}/{capture_id}/by_chrom/{chrom}.plp.gz",
        sample_map=rules.merge_cel.output.sample_map,
        snp_map=rules.merge_var.output.snp_map,
    output:
        plp=temp(out_dir / "popscle/{platform}/{capture_id}/wrangle/{chrom}.plp"),
    log:
        logs / "popscle/{platform}/{capture_id}/wrangle_plp_{chrom}.log",
    script:
        "../scripts/wrangle_plp.py"


rule merge_plp:
    input:
        plp=expand(
            out_dir / "popscle/{{platform}}/{{capture_id}}/wrangle/{chrom}.plp",
            chrom=chromosomes,
        ),
    output:
        plp=out_dir / "popscle/{platform}/{capture_id}/pileup.plp.gz",
    log:
        logs / "popscle/{platform}/{capture_id}/merge_plp.log",
    container:
        config["containers"]["bcftools"]
    shell:
        # Stack the per-chromosome wrangled PLPs (identical tab-separated schema)
        # keeping a single header. Replaces snp_demux's csvkit dependency (an
        # external docker image) with a header-dedup concat in the local
        # container — equivalent because every chunk shares the wrangle schema.
        r"""
        first=1
        for f in {input.plp}; do
            if [ "$first" -eq 1 ]; then
                cat "$f"
                first=0
            else tail -n +2 "$f"; fi
        done | gzip -n >{output.plp} 2>{log}
        """


rule wrangle_umi_by_chromosome:
    input:
        umi=out_dir / "popscle/{platform}/{capture_id}/by_chrom/{chrom}.umi.gz",
        sample_map=rules.merge_cel.output.sample_map,
    output:
        umi=temp(out_dir / "popscle/{platform}/{capture_id}/wrangle/{chrom}.umi.gz"),
    log:
        logs / "popscle/{platform}/{capture_id}/wrangle_umi_{chrom}.log",
    script:
        "../scripts/wrangle_umi.py"


rule merge_umi:
    input:
        umi=expand(
            out_dir / "popscle/{{platform}}/{{capture_id}}/wrangle/{chrom}.umi.gz",
            chrom=chromosomes,
        ),
    output:
        umi=out_dir / "popscle/{platform}/{capture_id}/pileup.umi.gz",
    log:
        logs / "popscle/{platform}/{capture_id}/merge_umi.log",
    shell:
        # -n: strip the gzip mtime so pileup.umi.gz is byte-reproducible (matches
        # merge_plp; keeps this stage's output stable in a DVC superpipeline).
        "zcat {input.umi} | gzip -n > {output.umi}"
