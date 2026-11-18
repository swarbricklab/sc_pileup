"""Summary QC metrics for one capture's pileup outputs (goal-1: additive,
reads the engines' existing outputs, changes nothing downstream).

Emits a tidy CSV: one row per engine x barcode set that ran for this capture
(barcode_set `cells`, plus a cellsnp `empties` row for soup pileups), with
coverage metrics. Deterministic (plain CSV, fixed column/row order, rounded floats).
"""
import csv
import gzip

import pandas as pd

snakemake = snakemake  # noqa: F821  (injected by Snakemake)


def parse_mtx(path):
    """Stream a MatrixMarket coordinate file. Returns dims + sums over nonzero
    entries. Handles both tab- (cellSNP) and space- (souporcell) separators and
    integer/real values."""
    n_rows = n_cols = 0
    total = 0.0
    rows_covered, cols_covered = set(), set()
    header_seen = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if not parts:
                continue
            if not header_seen:
                n_rows, n_cols = int(parts[0]), int(parts[1])
                header_seen = True
                continue
            r, c, v = int(parts[0]), int(parts[1]), float(parts[2])
            if v != 0:
                total += v
                rows_covered.add(r)
                cols_covered.add(c)
    return {
        "n_rows": n_rows, "n_cols": n_cols, "total": total,
        "n_rows_covered": len(rows_covered), "n_cols_covered": len(cols_covered),
    }


def gz_rows(path):
    """Number of data rows in a gzipped TSV (excludes the header line)."""
    with gzip.open(path, "rt") as fh:
        n = sum(1 for _ in fh)
    return max(n - 1, 0)


def inp(name):
    return getattr(snakemake.input, name, None)


capture_id = snakemake.wildcards.capture_id
platform = snakemake.wildcards.platform
rows = []

def cellsnp_row(prefix, barcode_set):
    """Metrics for one cellSNP-lite MEX: AD (alt), DP (total), OTH (other)
    sparse matrices (sites x cells)."""
    dp = parse_mtx(inp(prefix + "dp"))
    ad = parse_mtx(inp(prefix + "ad")) if inp(prefix + "ad") else {"total": None}
    oth = parse_mtx(inp(prefix + "oth")) if inp(prefix + "oth") else {"total": None}
    total_dp = dp["total"]
    total_alt = ad["total"]
    return {
        "capture_id": capture_id, "platform": platform, "engine": "cellsnp",
        "barcode_set": barcode_set,
        "n_cells": dp["n_cols"], "n_sites": dp["n_rows"],
        "n_sites_covered": dp["n_rows_covered"], "n_cells_covered": dp["n_cols_covered"],
        "total_depth": total_dp, "total_alt": total_alt, "total_other": oth["total"],
        "alt_fraction": (total_alt / total_dp) if total_dp else None,
        "mean_depth_per_covered_site": (total_dp / dp["n_rows_covered"]) if dp["n_rows_covered"] else None,
    }


if inp("cellsnp_dp"):
    rows.append(cellsnp_row("cellsnp_", "cells"))

# Soup pileup over the empty droplets (barcode_sets: empties; cellsnp only).
# n_cells here counts empty-droplet barcodes; umis_per_barcode_median comes from
# select_empty_barcodes' summary (Gene Expression UMIs per selected barcode).
if inp("empties_cellsnp_dp"):
    row = cellsnp_row("empties_cellsnp_", "empties")
    summ = pd.read_csv(inp("empties_summary")).iloc[0]
    row["umis_per_barcode_median"] = summ["umi_median"]
    rows.append(row)

# souporcell/vartrix: ref + alt sparse matrices (sites x cells)
if inp("souporcell_ref"):
    ref = parse_mtx(inp("souporcell_ref"))
    alt = parse_mtx(inp("souporcell_alt")) if inp("souporcell_alt") else {"total": 0.0, "n_rows_covered": 0}
    total_dp = ref["total"] + alt["total"]
    total_alt = alt["total"]
    n_cov = max(ref["n_rows_covered"], alt.get("n_rows_covered", 0))
    rows.append({
        "capture_id": capture_id, "platform": platform, "engine": "souporcell",
        "barcode_set": "cells",
        "n_cells": ref["n_cols"], "n_sites": ref["n_rows"],
        "n_sites_covered": n_cov, "n_cells_covered": ref["n_cols_covered"],
        "total_depth": total_dp, "total_alt": total_alt, "total_other": None,
        "alt_fraction": (total_alt / total_dp) if total_dp else None,
        "mean_depth_per_covered_site": (total_dp / n_cov) if n_cov else None,
    })

# popscle dsc-pileup: variant list + cell list (depth lives in plp, format-specific
# — reported in a later iteration). Counts only for now.
if inp("popscle_var"):
    rows.append({
        "capture_id": capture_id, "platform": platform, "engine": "popscle",
        "barcode_set": "cells",
        "n_cells": gz_rows(inp("popscle_cel")) if inp("popscle_cel") else None,
        "n_sites": gz_rows(inp("popscle_var")),
        "n_sites_covered": None, "n_cells_covered": None,
        "total_depth": None, "total_alt": None, "total_other": None,
        "alt_fraction": None, "mean_depth_per_covered_site": None,
    })

cols = ["capture_id", "platform", "engine", "barcode_set", "n_cells", "n_sites",
        "n_sites_covered", "n_cells_covered", "total_depth", "total_alt",
        "total_other", "alt_fraction", "mean_depth_per_covered_site"]
FLOAT_COLS = ["alt_fraction", "mean_depth_per_covered_site"]
# Soup-only metric: the column exists only when empties are piled up, so the
# default (cells-only) metrics differ from before by the barcode_set column alone.
if snakemake.params.empties:
    cols.append("umis_per_barcode_median")
    FLOAT_COLS.append("umis_per_barcode_median")
df = (pd.DataFrame(rows, columns=cols)
      .sort_values(["engine", "barcode_set"]).reset_index(drop=True))
# nullable Int64 so counts render as integers (not "3828.0") even when other
# engines leave the cell blank
INT_COLS = ["n_cells", "n_sites", "n_sites_covered", "n_cells_covered",
            "total_depth", "total_alt", "total_other"]
for c in INT_COLS:
    df[c] = df[c].round().astype("Int64")
# round floats for cross-platform-stable, churn-free CSV
for c in FLOAT_COLS:
    df[c] = pd.to_numeric(df[c]).round(6)
df.to_csv(snakemake.output[0], index=False, lineterminator="\n",
          quoting=csv.QUOTE_MINIMAL)
