#!/usr/bin/env python
"""Track-A validation: do cells carrying somatic ALT reads concentrate in the
epithelial/abnormal compartment? Joins the somatic cellSNP pileup to
cellranger's per-capture cell_types.csv (cheap; no h5ad)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from compare_pileups import load_cellsnp, _norm_chrom

RES = "/scratch/a56/jr9959/sc_pileup/results/cellsnp/short/Neomet_Pool028a_somatic"
MAP = "/scratch/a56/jr9959/sc_pileup/resources/somatic/site_donor_map.tsv"
CT = "/g/data/a56/dvc/registries/bcarc_chromium/files/md5/51/1f39564b245b032d723e1b84473ed8"  # cell_types.csv

# editing-prone site set
edit = set()
with open(MAP) as fh:
    next(fh)
    for line in fh:
        c, p, r, a, donors, ep = line.rstrip("\n").split("\t")
        if ep == "1":
            edit.add((_norm_chrom(c), p, r, a))

cs = load_cellsnp(RES)
sites, cells = cs["sites"], cs["cells"]
alt = cs["alt"].tocsr()
keep = np.array([s not in edit for s in sites])          # exclude RNA-editing-prone
burden = np.asarray((alt[keep, :] > 0).sum(axis=0)).ravel()   # somatic-alt sites per cell

# cell type per barcode
ctype = {}
with open(CT) as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split(",")
        ctype[f[0]] = f[1]
ct = np.array([ctype.get(b, "NA") for b in cells])

print(f"cells={len(cells)}  somatic-alt cells (>=1, excl editing)={int((burden>=1).sum())} "
      f"({(burden>=1).mean():.1%})\n")
print(f"{'coarse_cell_type':>22} {'n':>6} {'%>=1 alt':>9} {'%>=2':>7} {'mean burden':>12} {'enrich':>7}")
overall = (burden >= 1).mean()
rows = []
for t in sorted(set(ct)):
    m = ct == t
    if m.sum() < 20:
        continue
    p1 = (burden[m] >= 1).mean()
    rows.append((p1, t, m.sum(), (burden[m] >= 2).mean(), burden[m].mean()))
for p1, t, n, p2, mb in sorted(rows, reverse=True):
    print(f"{t:>22} {n:>6} {p1:>8.1%} {p2:>6.1%} {mb:>12.3f} {p1/overall:>6.1f}x")

# composition of the somatic-positive cells
print("\ncomposition of somatic-ALT+ cells (>=1, excl editing):")
pos = burden >= 1
for t, n in sorted([(t, (ct[pos] == t).sum()) for t in set(ct[pos])], key=lambda x: -x[1])[:8]:
    print(f"  {t:>22}: {n} ({n/pos.sum():.1%})")
