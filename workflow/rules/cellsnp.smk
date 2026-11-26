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
#   * per-platform extra flags via cellsnp_extra(), and per-platform
#     --minMAPQ/--minLEN via cellsnp_platform_param() (STAR vs minimap2 MAPQ);
#   * optional nM/NH read pre-filter (params.read_filter; rules/read_filter.smk).


wildcard_constraints:
    platform="short|long",


rule cellsnp_pileup:
    input:
        bam=cellsnp_bam,
        bai=cellsnp_bai,
        barcodes=lookup("barcodes"),
        vcf=sites_vcf_input,
        site_check=site_check_input,
    output:
        base=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.base.vcf.gz",
        samples=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.samples.tsv",
        ad=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.tag.AD.mtx",
        dp=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.tag.DP.mtx",
        oth=out_dir / "cellsnp/{platform}/{capture_id}/cellSNP.tag.OTH.mtx",
    log:
        logs / "cellsnp/{platform}/{capture_id}.log",
    container:
        config["containers"]["cellsnp"]
    params:
        minMAF=config["params"]["cellsnp_lite"]["minMAF"],
        minCOUNT=config["params"]["cellsnp_lite"]["minCOUNT"],
        min_mapq=cellsnp_platform_param("min_mapq"),
        min_len=cellsnp_platform_param("min_len"),
        extra=cellsnp_extra,
    shell:
        # --cellTAG CB, --UMItag UR and --maxDEPTH 0 are cellsnp-lite 1.2.2's own
        # defaults with a barcode list (--UMItag Auto means UR, the raw UMI),
        # passed explicitly to guard the UMI tag: cellsnp-lite 1.2.3 switches
        # that default to UB. --exclFLAG is left to cellsnp-lite, which derives
        # it from the UMI mode (duplicate-flagged reads are excluded only when no
        # UMI tag is used). Later releases also change how --maxDEPTH and
        # --maxPILEUP cap the pileup, so a version bump still needs a pileup
        # comparison at deep sites. platform_args ({params.extra}) come last on
        # the line, so a flag there, such as --UMItag UB, overrides these.
        """
        outdir=$(dirname {output.base})
        cellsnp-lite --samFile {input.bam} \
            --regionsVCF {input.vcf} \
            --barcodeFile {input.barcodes} \
            --outDir $outdir \
            --nproc {resources.threads} \
            --cellTAG CB \
            --UMItag UR \
            --maxDEPTH 0 \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --minMAPQ {params.min_mapq} \
            --minLEN {params.min_len} \
            --gzip {params.extra} \
            >{log} 2>&1
        """


# Same pileup over the capture's empty-droplet barcodes (barcode_sets: empties),
# for soup genotyping. Output under <outs.root>/empties/cellsnp/, mirroring the
# cells layout, so the cells paths are untouched and nothing nests inside them.
use rule cellsnp_pileup as cellsnp_pileup_empties with:
    input:
        bam=cellsnp_bam,
        bai=cellsnp_bai,
        barcodes=rules.select_empty_barcodes.output.barcodes,
        vcf=sites_vcf_input,
        site_check=site_check_input,
    output:
        base=out_dir / "empties/cellsnp/{platform}/{capture_id}/cellSNP.base.vcf.gz",
        samples=out_dir / "empties/cellsnp/{platform}/{capture_id}/cellSNP.samples.tsv",
        ad=out_dir / "empties/cellsnp/{platform}/{capture_id}/cellSNP.tag.AD.mtx",
        dp=out_dir / "empties/cellsnp/{platform}/{capture_id}/cellSNP.tag.DP.mtx",
        oth=out_dir / "empties/cellsnp/{platform}/{capture_id}/cellSNP.tag.OTH.mtx",
    log:
        logs / "empties/cellsnp/{platform}/{capture_id}.log",
