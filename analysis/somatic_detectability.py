#!/usr/bin/env python
"""Track-A feasibility: how many somatic variants are detectable in the pool's
scRNA, per donor, and what per-cell somatic-alt burden do we see?

Reads the somatic cellSNP pileup (minMAF 0/minCOUNT 1, so cellSNP only emits
sites with >=1 UMI => "detectable") + the site->donor map.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from compare_pileups import load_cellsnp, _norm_chrom

RES = "/scratch/a56/jr9959/sc_pileup/results/cellsnp/short/Neomet_Pool028a_somatic"
MAP = "/scratch/a56/jr9959/sc_pileup/resources/somatic/site_donor_map.tsv"

# site -> (donors set, editing_prone) ; per-donor input counts
donor_of = {}
edit_of = {}
input_per_donor = {}
n_input = 0
with open(MAP) as fh:
    next(fh)
    for line in fh:
        c, p, r, a, donors, ep = line.rstrip("\n").split("\t")
        k = (_norm_chrom(c), p, r, a)
        ds = donors.split(",")
        donor_of[k] = ds
        edit_of[k] = ep == "1"
        n_input += 1
        for d in ds:
            input_per_donor[d] = input_per_donor.get(d, 0) + 1

cs = load_cellsnp(RES)
sites, cells = cs["sites"], cs["cells"]
alt, total = cs["alt"].tocsr(), cs["total"].tocsr()          # sites x cells
site_alt = np.asarray(alt.sum(axis=1)).ravel()               # alt reads per site
n_det = len(sites)
print(f"somatic sites: input={n_input}  detected(>=1 UMI)={n_det} ({n_det/n_input:.1%})")

# classify detected sites
det_donor = np.array([donor_of.get(s, ["?"])[0] for s in sites])  # 1 donor each (private)
det_edit = np.array([edit_of.get(s, False) for s in sites])
has_alt = site_alt > 0

print("\n== per-donor detectability ==")
print(f"  {'donor':>6} {'input':>7} {'detected':>9} {'alt>0':>7} {'alt>0(noEdit)':>13}")
for d in sorted(input_per_donor):
    m = det_donor == d
    ae = m & has_alt & (~det_edit)
    print(f"  {d:>6} {input_per_donor[d]:>7} {m.sum():>9} {(m&has_alt).sum():>7} {ae.sum():>13}")

print("\n== alt-detected somatic sites (the usable Track-A signal) ==")
for lab, keep in [("all", np.ones(n_det, bool)), ("excl A>G/T>C editing", ~det_edit)]:
    k = keep & has_alt
    print(f"  {lab:>22}: {k.sum()} sites with >=1 somatic ALT read")

print("\n== per-cell somatic-ALT burden (candidate-malignant signal) ==")
for lab, keep in [("all sites", np.ones(n_det, bool)), ("excl editing", ~det_edit)]:
    sub = alt[keep, :]
    per_cell = np.asarray((sub > 0).sum(axis=0)).ravel()     # #somatic-alt sites per cell
    print(f"  {lab}: cells with >=1 somatic ALT: {(per_cell>=1).sum()} "
          f"({(per_cell>=1).mean():.1%}); >=2: {(per_cell>=2).sum()}; "
          f">=3: {(per_cell>=3).sum()}; max/cell: {per_cell.max()}")
print(f"\n  (n cells = {len(cells)})")
