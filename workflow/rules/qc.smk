# qc.smk — additive summary QC metrics for the pileup outputs.
#
# Goal-1 / downstream-neutral: these rules only READ the engines' existing
# outputs and emit CSV summaries; they add nothing to, and change nothing in,
# the pileup data the demux/somatic stages consume. Gated by `qc:` in config
# (default on). Cohort CSV: one row per capture x engine.


def _capture_metric_inputs(wildcards):
    """The engine outputs to summarise for this capture (only engines that ran)."""
    cid = wildcards.capture_id
    plat = captures.loc[cid, "platform"]
    files = {}
    if "cellsnp" in engines:
        base = out_dir / "cellsnp" / plat / cid
        files["cellsnp_dp"] = str(base / "cellSNP.tag.DP.mtx")
        files["cellsnp_ad"] = str(base / "cellSNP.tag.AD.mtx")
        files["cellsnp_oth"] = str(base / "cellSNP.tag.OTH.mtx")
    if "popscle" in engines and plat in _popscle_platforms:
        base = out_dir / "popscle" / plat / cid
        files["popscle_var"] = str(base / "pileup.var.gz")
        files["popscle_cel"] = str(base / "pileup.cel.gz")
    if "souporcell" in engines and plat in _souporcell_platforms:
        base = out_dir / "souporcell" / plat / cid
        files["souporcell_ref"] = str(base / "ref.mtx")
        files["souporcell_alt"] = str(base / "alt.mtx")
    return files


rule capture_metrics:
    input:
        unpack(_capture_metric_inputs),
    output:
        out_dir / "qc/{platform}/{capture_id}.metrics.csv",
    log:
        logs / "qc/{platform}/{capture_id}.metrics.log",
    script:
        "../scripts/capture_metrics.py"


rule cohort_metrics:
    input:
        lambda wc: [
            str(
                out_dir
                / "qc"
                / captures.loc[c, "platform"]
                / (str(c) + ".metrics.csv")
            )
            for c in captures.index
        ],
    output:
        out_dir / "qc/cohort_metrics.csv",
    log:
        logs / "qc/cohort_metrics.log",
    script:
        "../scripts/merge_metrics.py"
