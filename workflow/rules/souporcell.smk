# souporcell prep engine: minimap2 remap + vartrix -> ref.mtx / alt.mtx.
#
# Owns the expensive, deterministic, BAM-reading half of souporcell so snp_demux
# can resume from the matrices (clustering + troublet + consensus). See
# docs/souporcell_split.md.
#
# Genotype-free: de-novo variant calling with freebayes (souporcell's
# no---common_variants path), then vartrix counts at those variants. This matches
# the old snp_demux run, which ran souporcell_pipeline.py without --common_variants.
# The cluster count k and any cluster->donor mapping (--known_genotypes) stay in
# snp_demux.
#
# Commands are lifted verbatim from souporcell_pipeline.py (renamer -> minimap2
# -> retag -> freebayes -> vartrix); the whole engine runs in
# cumulusprod/souporcell:2.5, which carries minimap2, renamer/retag, samtools,
# bcftools, freebayes and vartrix.
#
# v1 is whole-BAM. minimap2 remap is not cleanly chromosome-splittable (reads
# relocate genome-wide); the depth+vartrix half is, and is deferred to v1.1.
#
# skip_remap (DEFAULT true): vartrix counts on the ORIGINAL cellranger BAM
# instead of the minimap2-remapped BAM. Validated to be a no-op: vartrix
# original-vs-remapped at somatic sites differs by ~1% in total ALT reads and
# 98.7% of per-cell malignant calls are unchanged, so the ~6 h + ~285 GB remap
# buys nothing. Set skip_remap: false to run the full remap. See
# docs/souporcell_split.md.

# Location of renamer.py / retag.py inside the souporcell image.
SOUPORCELL_DIR = "/opt/souporcell"


def _skip_remap():
    return config["params"]["souporcell"].get("skip_remap", True)


def souporcell_bam(wildcards):
    """BAM that vartrix counts on: original cellranger BAM (skip_remap) or the
    minimap2-remapped BAM (opt-in). When skip_remap, the remap rule is never
    requested and so never runs."""
    if _skip_remap():
        return captures.loc[wildcards.capture_id, "bam"]
    return str(
        out_dir
        / "souporcell"
        / wildcards.platform
        / wildcards.capture_id
        / "remap/minimap_tagged_sorted.bam"
    )


def souporcell_bai(wildcards):
    return souporcell_bam(wildcards) + ".bai"


rule souporcell_remap:
    input:
        bam=lookup("bam"),
        barcodes=lookup("barcodes"),
        fasta=config["ref"]["souporcell"]["genome_fasta"],
    output:
        # NOT temp(): the remapped BAM is a reusable asset (somatic Track-A
        # pileup on the remapped bam, re-vartrix, QC) and the remap is the
        # expensive step (~6 h). ~24 GB/capture — revisit a keep/discard config
        # flag if kept cohort-wide.
        bam=out_dir
        / "souporcell/{platform}/{capture_id}/remap/minimap_tagged_sorted.bam",
        bai=out_dir
        / "souporcell/{platform}/{capture_id}/remap/minimap_tagged_sorted.bam.bai",
    log:
        logs / "souporcell/{platform}/{capture_id}/remap.log",
    container:
        config["containers"]["souporcell"]
    params:
        soup=SOUPORCELL_DIR,
    shell:
        r"""
        d=$(dirname {output.bam})
        mkdir -p $d
        # 1) embed CB/UB into read names, emit fastq (whole BAM, single pass)
        {params.soup}/renamer.py --bam {input.bam} --barcodes {input.barcodes} \
            --out $d/reads.fq --no_umi False --umi_tag UB --cell_tag CB >{log} 2>&1
        # 2) splice-aware remap with souporcell's tuned minimap2 parameters
        minimap2 -ax splice -t {resources.threads} -G50k -k 21 -w 11 --sr \
            -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 \
            -g2000 -2K50m --secondary=no {input.fasta} $d/reads.fq \
            >$d/minimap.sam 2>>{log}
        # 3) restore CB/UB tags on the remapped reads
        {params.soup}/retag.py --sam $d/minimap.sam --no_umi False \
            --umi_tag UB --cell_tag CB --out $d/tagged.bam >>{log} 2>&1
        # 4) sort + index
        samtools sort -@ {resources.threads} -o {output.bam} $d/tagged.bam 2>>{log}
        samtools index {output.bam} 2>>{log}
        rm -f $d/reads.fq $d/minimap.sam $d/tagged.bam
        """


rule souporcell_call_variants_by_chromosome:
    input:
        bam=souporcell_bam,
        bai=souporcell_bai,
        fasta=config["ref"]["souporcell"]["genome_fasta"],
    output:
        vcf=temp(out_dir / "souporcell/{platform}/{capture_id}/by_chrom/{chrom}.vcf"),
    log:
        logs / "souporcell/{platform}/{capture_id}/call_variants_{chrom}.log",
    container:
        config["containers"]["souporcell"]
    params:
        min_coverage=int(config["params"]["souporcell"]["min_alt"])
        + int(config["params"]["souporcell"]["min_ref"]),
    shell:
        # De-novo variant calling, faithful to souporcell_pipeline.py's freebayes()
        # on the no---common_variants path: the exact flags souporcell runs (with
        # --min-coverage = min_alt + min_ref). Split per autosome (the `chromosomes`
        # list popscle already uses): souporcell only region-splits for parallelism,
        # so this is the same calling, just parallel and bounded per job. The alt /
        # decoy / sex contigs are skipped on purpose -- the imputed VCF carries only
        # chr1-22,X,Y, so a de-novo call off the autosomes can never match a donor
        # SNP downstream (and scRNA coverage there is negligible anyway).
        r"""
        freebayes -f {input.fasta} -iXu -C 2 -q 20 -n 3 -E 1 -m 30 \
            --min-coverage {params.min_coverage} --pooled-continuous --skip-coverage 100000 \
            -r {wildcards.chrom} {input.bam} >{output.vcf} 2>{log}
        """


rule souporcell_merge_variants:
    input:
        vcfs=expand(
            out_dir / "souporcell/{{platform}}/{{capture_id}}/by_chrom/{chrom}.vcf",
            chrom=chromosomes,
        ),
    output:
        vcf=temp(out_dir / "souporcell/{platform}/{capture_id}/variants_called.vcf"),
    log:
        logs / "souporcell/{platform}/{capture_id}/merge_variants.log",
    container:
        config["containers"]["souporcell"]
    shell:
        # Combine the per-chrom freebayes VCFs (identical headers) into one, then
        # sort to the reference contig order. Build header-once + records manually
        # (robust to the autosome-vs-fai header ordering) rather than bcftools concat.
        r"""
        tmp=$(dirname {output.vcf})/variants_unsorted.vcf
        bcftools view -h {input.vcfs[0]} >$tmp 2>{log}
        for f in {input.vcfs}; do bcftools view -H $f >>$tmp 2>>{log}; done
        bcftools sort -Ov -o {output.vcf} $tmp 2>>{log}
        rm -f $tmp
        """


rule souporcell_vartrix:
    input:
        bam=souporcell_bam,
        bai=souporcell_bai,
        vcf=rules.souporcell_merge_variants.output.vcf,
        barcodes=lookup("barcodes"),
        fasta=config["ref"]["souporcell"]["genome_fasta"],
    output:
        ref=out_dir / "souporcell/{platform}/{capture_id}/ref.mtx",
        alt=out_dir / "souporcell/{platform}/{capture_id}/alt.mtx",
        barcodes=out_dir / "souporcell/{platform}/{capture_id}/barcodes.tsv",
        variants=out_dir / "souporcell/{platform}/{capture_id}/variants.vcf",
    log:
        logs / "souporcell/{platform}/{capture_id}/vartrix.log",
    container:
        # vartrix (and the whole souporcell engine) runs in this image —
        # Demuxafy's bundled vartrix panics intermittently, which is why
        # snp_demux pins cumulusprod/souporcell:2.5.
        config["containers"]["souporcell"]
    params:
        mapq=config["params"]["souporcell"]["mapq"],
        extra=config["params"]["souporcell"]["extra"],
    shell:
        # Keep ref/alt matrices WITH their barcodes + variant VCF (the indices
        # that make them usable) — fixes snp_demux issue #35's orphaned matrices.
        r"""
        gunzip -c {input.barcodes} >{output.barcodes}
        cp {input.vcf} {output.variants}
        vartrix --mapq {params.mapq} -b {input.bam} -c {output.barcodes} \
            --scoring-method coverage --threads {resources.threads} \
            --ref-matrix {output.ref} --out-matrix {output.alt} \
            -v {output.variants} --fasta {input.fasta} --umi {params.extra} \
            >{log} 2>&1
        """
