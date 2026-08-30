# cellSNP-lite pileup — the primary, reusable artifact.
#
# Genotype-free: pile up every barcode at every site in sites_vcf. The MEX
# matrices (cell x SNP) feed vireo (demux) and, when sites_vcf is a somatic /
# CNV-region panel, Track A / Track B / cardelino.
#
# Reviewed vs snp_demux `rule cellsnp_lite`:
#   * outputs are NOT temp() (they are the product of this stage);
#   * site VCF is `sites_vcf` (generic), not the hard-coded donor panel;
#   * output path carries {platform} so short- and long-read pileups of the same
#     capture coexist and can be merged per cell downstream;
#   * per-platform extra flags via cellsnp_extra().

wildcard_constraints:
    platform="short|long",

rule cellsnp_pileup:
    input:
        bam=lookup("bam"),
        barcodes=lookup("barcodes"),
        vcf=config["deps"]["sites_vcf"],
    output:
        base=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.base.vcf.gz",
        samples=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.samples.tsv",
        ad=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.tag.AD.mtx",
        dp=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.tag.DP.mtx",
        oth=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.tag.OTH.mtx",
    params:
        minMAF=config["params"]["cellsnp_lite"]["minMAF"],
        minCOUNT=config["params"]["cellsnp_lite"]["minCOUNT"],
        extra=cellsnp_extra,
    threads: 32
    log:
        logs / "cellsnp/{platform}/{capture_id}.log",
    container:
        config["containers"]["demuxafy"]
    shell:
        """
        outdir=$(dirname {output.base})
        cellsnp-lite --samFile {input.bam} \
            --regionsVCF {input.vcf} \
            --barcodeFile {input.barcodes} \
            --outDir $outdir \
            --nproc {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --gzip {params.extra} \
            > {log} 2>&1
        """
