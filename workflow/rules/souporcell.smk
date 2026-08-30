# souporcell prep engine: minimap2 remap + vartrix -> ref.mtx / alt.mtx.
#
# Owns the expensive, deterministic, BAM-reading half of souporcell so snp_demux
# can resume from the matrices (clustering + troublet + consensus). See
# docs/souporcell_split.md.
#
# Genotype-free: vartrix counts at the covered subset of sites_vcf (souporcell's
# --common_variants path). The cluster count k and any cluster->donor mapping
# (--known_genotypes) stay in snp_demux.
#
# Commands are lifted verbatim from souporcell_pipeline.py (renamer -> minimap2
# -> retag -> covered-variant depth filter -> vartrix); pinned to the toolchain
# in Demuxafy.sif (minimap2 2.7, vartrix 1.1.3), whose prep is equivalent to
# cumulusprod/souporcell:2.5.
#
# v1 is whole-BAM. minimap2 remap is not cleanly chromosome-splittable (reads
# relocate genome-wide); the depth+vartrix half is, and is deferred to v1.1.

# Location of renamer.py / retag.py inside Demuxafy.sif.
SOUPORCELL_DIR = "/opt/souporcell"


rule souporcell_remap:
    input:
        bam=lookup("bam"),
        barcodes=lookup("barcodes"),
        fasta=config["ref"]["souporcell"]["genome_fasta"],
    output:
        bam=temp(out_dir / "souporcell/{platform}/{capture_id}/remap/minimap_tagged_sorted.bam"),
        bai=temp(out_dir / "souporcell/{platform}/{capture_id}/remap/minimap_tagged_sorted.bam.bai"),
    params:
        soup=SOUPORCELL_DIR,
    log:
        logs / "souporcell/{platform}/{capture_id}/remap.log",
    container:
        config["containers"]["demuxafy"]
    shell:
        r"""
        d=$(dirname {output.bam}); mkdir -p $d
        # 1) embed CB/UB into read names, emit fastq (whole BAM, single pass)
        {params.soup}/renamer.py --bam {input.bam} --barcodes {input.barcodes} \
            --out $d/reads.fq --no_umi False --umi_tag UB --cell_tag CB > {log} 2>&1
        # 2) splice-aware remap with souporcell's tuned minimap2 parameters
        minimap2 -ax splice -t {resources.threads} -G50k -k 21 -w 11 --sr \
            -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 \
            -g2000 -2K50m --secondary=no {input.fasta} $d/reads.fq \
            > $d/minimap.sam 2>> {log}
        # 3) restore CB/UB tags on the remapped reads
        {params.soup}/retag.py --sam $d/minimap.sam --no_umi False \
            --umi_tag UB --cell_tag CB --out $d/tagged.bam >> {log} 2>&1
        # 4) sort + index
        samtools sort -@ {resources.threads} -o {output.bam} $d/tagged.bam 2>> {log}
        samtools index {output.bam} 2>> {log}
        rm -f $d/reads.fq $d/minimap.sam $d/tagged.bam
        """


rule souporcell_covered_variants:
    input:
        bam=rules.souporcell_remap.output.bam,
        bai=rules.souporcell_remap.output.bai,
        vcf=config["deps"]["sites_vcf"],
    output:
        vcf=temp(out_dir / "souporcell/{platform}/{capture_id}/covered_variants.vcf"),
    params:
        min_cov=config["params"]["souporcell"]["min_cov"],
    log:
        logs / "souporcell/{platform}/{capture_id}/covered_variants.log",
    container:
        config["containers"]["demuxafy"]
    shell:
        # Keep only sites with read depth >= min_cov on the remapped bam, then
        # subset sites_vcf to those covered regions (mirrors souporcell). awk
        # form matches souporcell's ($2, $2+1) bed convention.
        r"""
        d=$(dirname {output.vcf})
        samtools depth {input.bam} \
          | awk -v c={params.min_cov} '{{ if ($3 >= c && $3 < 100000) print $1"\t"$2"\t"($2+1)"\t"$3 }}' \
          | bedtools merge -i - > $d/depth_merged.bed 2> {log}
        ( bcftools view -h {input.vcf}; \
          bedtools intersect -wa -a {input.vcf} -b $d/depth_merged.bed ) \
          > {output.vcf} 2>> {log}
        rm -f $d/depth_merged.bed
        """


rule souporcell_vartrix:
    input:
        bam=rules.souporcell_remap.output.bam,
        bai=rules.souporcell_remap.output.bai,
        vcf=rules.souporcell_covered_variants.output.vcf,
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
        # NB: cumulusprod/souporcell:2.5, NOT Demuxafy — Demuxafy's bundled
        # vartrix panics intermittently (this is why snp_demux pins souporcell
        # to the same image). remap/covered use the local Demuxafy tools (fine).
        config["containers"]["souporcell"]
    shell:
        # Keep ref/alt matrices WITH their barcodes + variant VCF (the indices
        # that make them usable) — fixes snp_demux issue #35's orphaned matrices.
        r"""
        gunzip -c {input.barcodes} > {output.barcodes}
        cp {input.vcf} {output.variants}
        vartrix --mapq 30 -b {input.bam} -c {output.barcodes} \
            --scoring-method coverage --threads {resources.threads} \
            --ref-matrix {output.ref} --out-matrix {output.alt} \
            -v {output.variants} --fasta {input.fasta} --umi \
            > {log} 2>&1
        """
