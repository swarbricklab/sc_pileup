# Optional per-capture BAM read pre-filter for the cellsnp engine.
#
# cellSNP-lite admits reads on MAPQ and mapped length only (--minMAPQ/--minLEN).
# For somatic sites, where the signal is often 1-2 UMIs per cell, multi-mismatch
# and multi-mapping reads are a meaningful share of false ALT counts, so
# `params.read_filter: {max_nM, max_NH}` (default off) inserts a
# `samtools view -e` pass whose output cellsnp reads in place of the original BAM
# (SComatic recommends max_nM 5, max_NH 1 for scRNA). Reads that lack the tag
# (e.g. minimap2 long-read BAMs carry no nM/NH) are KEPT -- the filter only
# removes reads that positively fail it. The log counts the reads lacking each
# capped tag, and warns when every read lacks it (that cap is then inert).
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
        tags=" ".join(read_filter_tags),
    shell:
        """
        samtools view -b -@ {resources.threads} -e '{params.expr}' \
            -o {output.bam} {input.bam} 2>{log}
        samtools index {output.bam} 2>>{log}
        echo "read filter: {params.expr}" >>{log}
        n_in=$(samtools view -c -@ {resources.threads} {input.bam})
        echo "reads in:  $n_in" >>{log}
        echo "reads out: $(samtools view -c -@ {resources.threads} {output.bam})" >>{log}
        for tag in {params.tags}; do
            n_lack=$(samtools view -c -@ {resources.threads} -e "![$tag]" {input.bam})
            echo "reads lacking $tag (kept): $n_lack" >>{log}
            if [ "$n_in" -gt 0 ] && [ "$n_lack" -eq "$n_in" ]; then
                echo "WARNING: no read in {input.bam} carries a $tag tag, so the" \
                    "$tag cap removes nothing" | tee -a {log} >&2
            fi
        done
        """
