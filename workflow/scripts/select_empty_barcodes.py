"""Select empty-droplet barcodes for one capture (soup pileups).

Empty droplets = barcodes in the Cell Ranger RAW matrix that were NOT called as
cells and whose total UMI count falls inside params.empties.umi_window
[min, max] (inclusive). Only "Gene Expression" features are counted when the
features file has a feature-type column (antibody / CRISPR capture UMIs would
otherwise inflate the totals).

Outputs:
  barcodes.tsv.gz  one barcode per line, sorted (the format cellsnp-lite takes)
  summary.csv      n_selected, UMI quantiles of the selected set, window used

Deterministic: sorted barcodes, gzip mtime fixed at 0, rounded floats.
"""
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

snakemake = snakemake  # noqa: F821  (injected by Snakemake)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

GEX = "Gene Expression"
CHUNK = 10_000_000


def first_existing(d, names):
    for n in names:
        if (d / n).exists():
            return d / n
    raise FileNotFoundError(f"none of {names} in {d}")


def read_lines(path):
    return pd.read_csv(path, sep="\t", header=None, dtype=str)


raw = Path(snakemake.input.raw_matrix)
lo, hi = (int(x) for x in snakemake.params.umi_window)
capture_id = snakemake.wildcards.capture_id

mtx = first_existing(raw, ["matrix.mtx.gz", "matrix.mtx"])
bcs = read_lines(first_existing(raw, ["barcodes.tsv.gz", "barcodes.tsv"]))[0].to_numpy()
feats = read_lines(first_existing(raw, ["features.tsv.gz", "features.tsv", "genes.tsv"]))
if feats.shape[1] >= 3:
    keep_feature = (feats[2] == GEX).to_numpy()
else:  # old-style genes.tsv: all gene expression
    keep_feature = np.ones(len(feats), dtype=bool)

# Stream the coordinate entries (1-based feature, barcode, count), summing
# per-barcode UMIs over the kept features.
with (gzip.open(mtx, "rt") if mtx.suffix == ".gz" else open(mtx)) as fh:
    n_header = 0
    for line in fh:
        n_header += 1
        if not line.startswith("%"):
            n_feat, n_bc, _ = (int(x) for x in line.split())
            break
if n_feat != len(feats) or n_bc != len(bcs):
    raise ValueError(
        f"{mtx}: dims {n_feat}x{n_bc} don't match {len(feats)} features x {len(bcs)} barcodes"
    )
totals = np.zeros(n_bc, dtype=np.int64)
for chunk in pd.read_csv(mtx, sep=r"\s+", header=None, skiprows=n_header,
                         dtype=np.int64, chunksize=CHUNK, engine="c"):
    f, b, v = (chunk[i].to_numpy() for i in (0, 1, 2))
    m = keep_feature[f - 1]
    totals += np.bincount(b[m] - 1, weights=v[m], minlength=n_bc).astype(np.int64)

cells = set(read_lines(snakemake.input.cell_barcodes)[0])
is_cell = np.fromiter((b in cells for b in bcs), dtype=bool, count=len(bcs))
in_window = (totals >= lo) & (totals <= hi)
sel = in_window & ~is_cell

selected = sorted(bcs[sel])
sel_umis = totals[sel]
with open(snakemake.output.barcodes, "wb") as raw_fh:
    with gzip.GzipFile(fileobj=raw_fh, mode="wb", mtime=0, filename="") as gz:
        gz.write("".join(b + "\n" for b in selected).encode())

q = (np.quantile(sel_umis, [0.05, 0.25, 0.5, 0.75, 0.95]) if len(sel_umis)
     else [None] * 5)
summary = pd.DataFrame([{
    "capture_id": capture_id,
    "umi_window_min": lo,
    "umi_window_max": hi,
    "n_raw_barcodes": int(n_bc),
    "n_cell_barcodes": len(cells),
    "n_cell_barcodes_in_window": int((in_window & is_cell).sum()),
    "n_selected": len(selected),
    "total_umis_selected": int(sel_umis.sum()),
    "umi_q05": q[0], "umi_q25": q[1], "umi_median": q[2],
    "umi_q75": q[3], "umi_q95": q[4],
}])
for c in ["umi_q05", "umi_q25", "umi_median", "umi_q75", "umi_q95"]:
    summary[c] = pd.to_numeric(summary[c]).round(3)
summary.to_csv(snakemake.output.summary, index=False, lineterminator="\n",
               quoting=csv.QUOTE_MINIMAL)
if not selected:
    print(f"WARNING: no empty-droplet barcodes selected for {capture_id} "
          f"in UMI window [{lo}, {hi}]")
print(summary.to_string(index=False))
