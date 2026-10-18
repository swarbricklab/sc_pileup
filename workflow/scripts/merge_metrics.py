"""Concatenate per-capture metric CSVs into one cohort summary CSV.

Deterministic: sorted by (capture_id, engine), fixed newline, no index.
"""
import csv

import pandas as pd

snakemake = snakemake  # noqa: F821

# read_csv re-infers dtypes (losing the per-capture Int64), so re-cast counts
INT_COLS = ["n_cells", "n_sites", "n_sites_covered", "n_cells_covered",
            "total_depth", "total_alt", "total_other"]
frames = [pd.read_csv(f) for f in snakemake.input]
cohort = pd.concat(frames, ignore_index=True)
cohort = cohort.sort_values(["capture_id", "engine"]).reset_index(drop=True)
for c in INT_COLS:
    cohort[c] = cohort[c].round().astype("Int64")
cohort.to_csv(snakemake.output[0], index=False, lineterminator="\n",
              quoting=csv.QUOTE_MINIMAL)
