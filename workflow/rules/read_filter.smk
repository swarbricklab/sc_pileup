# Optional per-capture BAM read pre-filter for the cellsnp engine.
#
# cellSNP-lite admits reads on MAPQ and mapped length only (--minMAPQ/--minLEN).
# For somatic sites, where the signal is often 1-2 UMIs per cell, multi-mismatch
# and multi-mapping reads are a meaningful share of false ALT counts, so
# `params.read_filter: {max_nM, max_NH}` (default off) inserts a
# `samtools view -e` pass whose output cellsnp reads in place of the original BAM
# (SComatic recommends max_nM 5, max_NH 1 for scRNA). Reads that lack the tag
# (e.g. minimap2 long-read BAMs carry no nM/NH) are KEPT -- the filter only
# removes reads that positively fail it.
#
# Off (both null / block absent): this rule is never requested and the cellsnp
# DAG is unchanged. popscle and souporcell always read the original BAM.
#
# Not here: base-quality masking. cellSNP-lite has no BQ floor; its effect is
# measured downstream by per-site QC (somatifind), not filtered here.


rule read_filter_bam:
    input:
        bam=lookup("bam"),
    output:
        bam=temp(out_dir / "read_filter/{platform}/{capture_id}/filtered.bam"),
        bai=temp(out_dir / "read_filter/{platform}/{capture_id}/filtered.bam.bai"),
    log:
        logs / "read_filter/{platform}/{capture_id}.log",
    container:
        config["containers"]["samtools"]
    params:
        expr=read_filter_expr,
    shell:
        """
        samtools view -b -@ {resources.threads} -e '{params.expr}' \
            -o {output.bam} {input.bam} 2>{log}
        samtools index {output.bam} 2>>{log}
        echo "read filter: {params.expr}" >>{log}
        echo "reads in:  $(samtools view -c {input.bam})" >>{log}
        echo "reads out: $(samtools view -c {output.bam})" >>{log}
        """
