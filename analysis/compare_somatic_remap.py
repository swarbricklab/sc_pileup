#!/usr/bin/env python
"""Does remap change SOMATIC detection? cellSNP at the somatic panel on the
original cellranger BAM vs the souporcell-remapped BAM (same tool, same sites).
Decides whether skip-remap is safe for the higher-stakes somatic case."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from compare_pileups import load_cellsnp

ORIG = "/scratch/a56/jr9959/sc_pileup/results/cellsnp/short/Neomet_Pool028a_somatic"
REMAP = "/scratch/a56/jr9959/sc_pileup/results/cellsnp/short/Neomet_Pool028a_somatic_remap"

A = load_cellsnp(ORIG)
B = load_cellsnp(REMAP)


def keymap(e):
    return {k: i for i, k in enumerate(e["sites"])}, {c: j for j, c in enumerate(e["cells"])}


ak, ac = keymap(A)
bk, bc = keymap(B)
common = [k for k in A["sites"] if k in bk]
print(f"somatic sites detected: original={len(A['sites'])}  remapped={len(B['sites'])}  shared={len(common)}")

# per-site alt totals
def site_alt(e):
    return np.asarray(e["alt"].tocsr().sum(axis=1)).ravel()

sa, sb = site_alt(A), site_alt(B)
oa = np.array([sa[ak[k]] for k in common]); ob = np.array([sb[bk[k]] for k in common])
alt_o = (oa > 0).sum(); alt_b = (ob > 0).sum()
print(f"somatic sites with >=1 ALT read: original={ (sa>0).sum() }  remapped={ (sb>0).sum() }  "
      f"(shared sites: orig {alt_o}, remap {alt_b})")
print(f"total somatic ALT reads (shared sites): original={int(oa.sum())}  remapped={int(ob.sum())}  "
      f"net change={int(ob.sum()-oa.sum()):+d} ({100*(ob.sum()-oa.sum())/max(oa.sum(),1):+.1f}%)")

# per-cell somatic burden (alt sites/cell) on shared cells
cells = [c for c in A["cells"] if c in bc]
ai = [ac[c] for c in cells]; bi = [bc[c] for c in cells]
si = [ak[k] for k in common]; sj = [bk[k] for k in common]
Asub = A["alt"].tocsr()[si, :][:, ai]
Bsub = B["alt"].tocsr()[sj, :][:, bi]
burden_o = np.asarray((Asub > 0).sum(axis=0)).ravel()
burden_b = np.asarray((Bsub > 0).sum(axis=0)).ravel()
for thr in (1, 2, 3):
    print(f"cells with >= {thr} somatic ALT: original={int((burden_o>=thr).sum())}  "
          f"remapped={int((burden_b>=thr).sum())}")
# do per-cell malignant calls (>=2) flip?
mo, mb = burden_o >= 2, burden_b >= 2
print(f"\nmalignant(>=2) cells: original={int(mo.sum())} remapped={int(mb.sum())}  "
      f"agree={ (mo==mb).mean():.1%}  gained={int((mb&~mo).sum())} lost={int((mo&~mb).sum())}")
