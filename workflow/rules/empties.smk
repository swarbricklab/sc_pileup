# Empty-droplet ("soup") barcode selection, per capture.
#
# Genotyping the ambient soup gives each capture's soup donor composition, used
# to scale genotype-detected cross-donor contamination up to total contamination
# and to ground SoupX/CellBender parameters. This picks which barcodes count as
# empty droplets: barcodes in the Cell Ranger RAW matrix (captures.csv
# `raw_matrix` column: the raw_feature_bc_matrix/ MEX directory) that were not
# called as cells (`barcodes` column) and whose Gene Expression UMI total is in
# params.empties.umi_window. See config.yaml for the window rationale.
#
# Only runs when its outputs are requested (barcode_sets: [.., empties]).


def raw_matrix_input(wildcards):
    if "raw_matrix" not in captures.columns:
        raise ValueError(
            "empty-droplet barcodes requested but captures.csv has no 'raw_matrix' "
            "column (path to the Cell Ranger raw_feature_bc_matrix/ directory)"
        )
    v = captures.loc[wildcards.capture_id, "raw_matrix"]
    if pd.isna(v) or not str(v).strip():
        raise ValueError(
            "empty-droplet barcodes requested but capture '"
            + str(wildcards.capture_id)
            + "' has an empty 'raw_matrix' value"
        )
    return str(v).strip()


rule select_empty_barcodes:
    input:
        raw_matrix=raw_matrix_input,
        cell_barcodes=lookup("barcodes"),
    output:
        barcodes=out_dir / "empties/{capture_id}/barcodes.tsv.gz",
        summary=out_dir / "empties/{capture_id}/summary.csv",
    log:
        logs / "empties/{capture_id}.log",
    params:
        umi_window=config["params"].get("empties", {}).get("umi_window", [10, 100]),
    script:
        "../scripts/select_empty_barcodes.py"
