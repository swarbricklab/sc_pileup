#!/usr/bin/env python
"""Characterise the cellSNP-vs-souporcell allele-fraction discordance:
is the ~3% tail low-coverage sampling noise, or systematic ref bias?

Decision logic:
  - if discordance collapses as the coverage floor rises  -> sampling noise -> remap unneeded
  - if it persists and is site-concentrated + directional -> real ref bias  -> remap (or site exclusion) matters
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from compare_pileups import load_cellsnp, load_souporcell, align

RES = "/scratch/a56/jr9959/sc_pileup/results"
CAP, PLAT = "Neomet_Pool028a", "short"

A = load_cellsnp(f"{RES}/cellsnp/{PLAT}/{CAP}")
B = load_souporcell(f"{RES}/souporcell/{PLAT}/{CAP}")
sites, cells, aA, tA, aB, tB = align(A, B)
aA, tA, aB, tB = aA.tocsr(), tA.tocsr(), aB.tocsr(), tB.tocsr()

# entries covered (>=5) in BOTH
both = (tA >= 5).multiply(tB >= 5).tocoo()
r, c = both.row, both.col
totA = np.asarray(tA[r, c]).ravel(); totB = np.asarray(tB[r, c]).ravel()
afA = np.asarray(aA[r, c]).ravel() / totA
afB = np.asarray(aB[r, c]).ravel() / totB
diff = afB - afA                      # souporcell(remapped) - cellsnp(original)
adiff = np.abs(diff)
disc = adiff > 0.1
covmin = np.minimum(totA, totB)
n = len(r)
print(f"shared covered(>=5) entries: {n};  discordant(|dAF|>0.1): {disc.sum()} ({disc.mean():.2%})")

print("\n== 1. coverage dependence (is it sampling noise?) ==")
for lo, hi in [(5,10),(10,20),(20,50),(50,10**9)]:
    m = (covmin >= lo) & (covmin < hi)
    lab = f"{lo}-{hi-1}" if hi < 10**9 else f"{lo}+"
    if m.sum():
        print(f"  min-cov {lab:>6}: n={m.sum():>7}  discordant={disc[m].mean():.2%}")

print("\n== 2. direction (ref bias suppresses ALT -> remap should RAISE alt frac) ==")
print(f"  mean signed dAF (all entries)       = {diff.mean():+.4f}")
print(f"  mean signed dAF (discordant only)   = {diff[disc].mean():+.4f}")
print(f"  discordant with souporcell>cellsnp  = {(diff[disc] > 0).mean():.1%}")

print("\n== 3. site concentration (few bad sites vs spread thin?) ==")
nsite = len(sites)
n_per_site = np.bincount(r, minlength=nsite).astype(float)
d_per_site = np.bincount(r, weights=disc.astype(float), minlength=nsite)
have = n_per_site >= 20                          # sites with >=20 covered cells
rate = np.divide(d_per_site, n_per_site, out=np.zeros(nsite), where=n_per_site>0)
print(f"  sites with >=20 covered cells: {have.sum()}")
for thr in [0.5, 0.25, 0.1]:
    hot = have & (rate >= thr)
    frac_disc = d_per_site[hot].sum() / disc.sum()
    print(f"  sites >= {thr:.0%} of cells discordant: {hot.sum():>5}  "
          f"(hold {frac_disc:.1%} of ALL discordant entries)")
# how concentrated: cumulative share from most-discordant sites
order = np.argsort(-d_per_site)
cum = np.cumsum(d_per_site[order]) / disc.sum()
for k in [50, 200, 1000]:
    print(f"  top {k:>4} sites by #discordant hold {cum[min(k,len(cum))-1]:.1%} of discordant entries")

print("\n== 4. are discordant sites the heterozygous ones (AF~0.5)? ==")
# per-site mean AF (cellsnp), weighted by covered cells
site_meanAF = np.zeros(nsite)
np.add.at(site_meanAF, r, afA)
site_meanAF = np.divide(site_meanAF, n_per_site, out=np.zeros(nsite), where=n_per_site>0)
hot = have & (rate >= 0.25)
def band(mask):
    v = site_meanAF[mask]
    return f"AF<0.15:{np.mean(v<0.15):.0%}  0.15-0.85:{np.mean((v>=0.15)&(v<=0.85)):.0%}  >0.85:{np.mean(v>0.85):.0%}"
print(f"  all covered sites : {band(have)}")
print(f"  hot(>=25%) sites  : {band(hot)}")

print("\n== 5. top 10 most-discordant sites (>=20 cells) ==")
top = [i for i in order if have[i]][:10]
print("  chrom:pos           n_cells  disc_rate  meanAF_cellsnp  meanAF_soup(mean signed)")
site_meanAFb = np.zeros(nsite); np.add.at(site_meanAFb, r, afB)
site_meanAFb = np.divide(site_meanAFb, n_per_site, out=np.zeros(nsite), where=n_per_site>0)
for i in top:
    ch, po, rf, al = sites[i]
    print(f"  {ch}:{po:<12} {int(n_per_site[i]):>7}  {rate[i]:>8.1%}  "
          f"{site_meanAF[i]:>13.3f}  {site_meanAFb[i]:>6.3f}")
