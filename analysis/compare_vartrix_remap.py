#!/usr/bin/env python
"""vartrix original(skip_remap) vs remapped, at somatic sites — the clean
skip-remap test. Same VCF+barcodes -> matrices align element-wise (no keying).
Editing-prone (A>G/T>C) sites excluded via the panel's site_donor_map."""
import gzip
import numpy as np
from scipy.io import mmread

OUT = "/scratch/a56/jr9959/sc_pileup/results/vartrix_compare/Neomet_Pool028a"
MAP = "/scratch/a56/jr9959/sc_pileup/resources/somatic/site_donor_map.tsv"
VCF = f"{OUT}/somatic.vcf"

# editing flag per site, in VCF row order
edit = {}
with open(MAP) as fh:
    next(fh)
    for l in fh:
        c, p, r, a, d, ep = l.rstrip("\n").split("\t")
        edit[(c.replace("chr", ""), p, r, a)] = ep == "1"
rows = []
with open(VCF) as fh:
    for l in fh:
        if l[0] == "#":
            continue
        f = l.split("\t")
        rows.append(edit.get((f[0].replace("chr", ""), f[1], f[3], f[4]), False))
keep = ~np.array(rows)  # non-editing

ao = mmread(f"{OUT}/alt_orig.mtx").tocsr()
ar = mmread(f"{OUT}/alt_remap.mtx").tocsr()
print(f"matrices: orig {ao.shape} nnz={ao.nnz}   remap {ar.shape} nnz={ar.nnz}")
ao, ar = ao[keep, :], ar[keep, :]

so = np.asarray(ao.sum(1)).ravel(); sr = np.asarray(ar.sum(1)).ravel()      # per-site alt
print(f"\nsomatic sites with >=1 ALT read:  orig={int((so>0).sum())}  remap={int((sr>0).sum())}")
print(f"total somatic ALT reads:          orig={int(so.sum())}  remap={int(sr.sum())}  "
      f"net={int(sr.sum()-so.sum()):+d} ({100*(sr.sum()-so.sum())/max(so.sum(),1):+.1f}%)")
print(f"sites gained by remap: {int(((so==0)&(sr>0)).sum())}   lost: {int(((so>0)&(sr==0)).sum())}")

bo = np.asarray((ao > 0).sum(0)).ravel(); br = np.asarray((ar > 0).sum(0)).ravel()  # per-cell burden
for t in (1, 2, 3):
    print(f"cells with >= {t} somatic ALT:  orig={int((bo>=t).sum())}  remap={int((br>=t).sum())}")
mo, mr = bo >= 2, br >= 2
print(f"\nmalignant(>=2): orig={int(mo.sum())} remap={int(mr.sum())}  agree={(mo==mr).mean():.1%}  "
      f"gained={int((mr&~mo).sum())} lost={int((mo&~mr).sum())}")
