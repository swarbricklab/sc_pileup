#!/usr/bin/env python
"""Cross-compare the sc_pileup engines (cellsnp / souporcell / popscle) for one
capture, putting them on a common (site x cell) grid and reporting agreement.

Two questions this answers:
  1. QC: do the engines agree on per-cell coverage and allele fractions at
     shared sites? (divergence = bug / format issue)
  2. Decision: cellSNP counts on the ORIGINAL bam, souporcell/vartrix on the
     REMAPPED bam — comparing their allele fractions quantifies souporcell's
     ref-bias correction, i.e. whether the separate remap engine earns its cost.

Engines are loaded best-effort: whatever is present under
  <results>/<engine>/<platform>/<capture>/
gets compared; missing engines are skipped (so this runs before souporcell
finishes, then again with all three).

cellSNP and souporcell give per-cell (alt, total) matrices -> full allele-
fraction comparison. popscle's plp is read-level; here it contributes only
coverage-level metrics (sites/UMIs per cell from cel.gz).

Usage:
  compare_pileups.py --results results --capture Neomet_Pool028a \
      --platform short --out results/compare/Neomet_Pool028a
"""
import argparse, gzip, os, sys, json
import numpy as np
import pandas as pd
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------- site-key helpers ----------
def _norm_chrom(c):
    return c[3:] if c.startswith("chr") else c


def _read_vcf_sites(path):
    """Return list of site keys (chrom,pos,ref,alt) in file (=matrix row) order."""
    op = gzip.open if path.endswith(".gz") else open
    sites = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            sites.append((_norm_chrom(f[0]), f[1], f[3], f[4]))
    return sites


def _read_lines(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return [ln.strip().split("\t")[0] for ln in fh if ln.strip()]


# ---------- per-engine loaders ----------
def load_cellsnp(d):
    sites = _read_vcf_sites(os.path.join(d, "cellSNP.base.vcf.gz"))
    cells = _read_lines(os.path.join(d, "cellSNP.samples.tsv"))
    alt = mmread(os.path.join(d, "cellSNP.tag.AD.mtx")).tocsr()   # sites x cells
    total = mmread(os.path.join(d, "cellSNP.tag.DP.mtx")).tocsr()
    return dict(name="cellsnp", sites=sites, cells=cells, alt=alt, total=total)


def load_souporcell(d):
    sites = _read_vcf_sites(os.path.join(d, "variants.vcf"))
    cells = _read_lines(os.path.join(d, "barcodes.tsv"))
    alt = mmread(os.path.join(d, "alt.mtx")).tocsr()              # variants x cells
    ref = mmread(os.path.join(d, "ref.mtx")).tocsr()
    return dict(name="souporcell", sites=sites, cells=cells, alt=alt, total=(alt + ref).tocsr())


def load_popscle(d):
    """Coverage-level only: per-cell sites/UMIs from cel.gz, site count from var.gz."""
    cel = pd.read_csv(os.path.join(d, "pileup.cel.gz"), sep="\t")
    nsite = sum(1 for _ in gzip.open(os.path.join(d, "pileup.var.gz"), "rt")) - 1
    cov = pd.Series(cel["NUM.SNP"].values, index=cel["BARCODE"].values)
    umi = pd.Series(cel["NUM.UMIwSNP"].values, index=cel["BARCODE"].values)
    return dict(name="popscle", n_sites=int(nsite), cells=list(cel["BARCODE"]),
                sites_per_cell=cov, umi_per_cell=umi)


# ---------- comparison of two matrix engines ----------
def align(a, b):
    """Common sites & cells between two matrix-engine dicts; return aligned
    dense-friendly CSR submatrices in identical (site x cell) order."""
    a_site = {k: i for i, k in enumerate(a["sites"])}
    b_site = {k: i for i, k in enumerate(b["sites"])}
    common_sites = [k for k in a["sites"] if k in b_site]
    a_cell = {k: i for i, k in enumerate(a["cells"])}
    b_cell = {k: i for i, k in enumerate(b["cells"])}
    common_cells = [k for k in a["cells"] if k in b_cell]
    ai = [a_site[k] for k in common_sites]; bi = [b_site[k] for k in common_sites]
    ac = [a_cell[k] for k in common_cells]; bc = [b_cell[k] for k in common_cells]
    sub = lambda m, r, c: m[r, :][:, c]
    return (common_sites, common_cells,
            sub(a["alt"], ai, ac), sub(a["total"], ai, ac),
            sub(b["alt"], bi, bc), sub(b["total"], bi, bc))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--capture", required=True)
    ap.add_argument("--platform", default="short")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cov", type=int, default=5,
                    help="min total count in BOTH engines for an AF comparison entry")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    root = lambda eng: os.path.join(args.results, eng, args.platform, args.capture)
    engines = {}
    loaders = dict(cellsnp=load_cellsnp, souporcell=load_souporcell, popscle=load_popscle)
    for eng, fn in loaders.items():
        d = root(eng)
        if os.path.isdir(d):
            try:
                engines[eng] = fn(d)
                print(f"loaded {eng} from {d}", file=sys.stderr)
            except Exception as e:
                print(f"skip {eng}: {e}", file=sys.stderr)

    report = {"capture": args.capture, "platform": args.platform, "engines": {}}
    for eng, e in engines.items():
        if "alt" in e:
            report["engines"][eng] = dict(n_sites=len(e["sites"]), n_cells=len(e["cells"]),
                                          nnz=int(e["total"].nnz))
        else:
            report["engines"][eng] = dict(n_sites=e["n_sites"], n_cells=len(e["cells"]))

    # ---- sites/cell per engine (coverage yield) ----
    fig, ax = plt.subplots(figsize=(6, 4))
    for eng, e in engines.items():
        if "total" in e:
            spc = np.asarray((e["total"] > 0).sum(axis=0)).ravel()
        else:
            spc = e["sites_per_cell"].values
        spc = spc[spc > 0]
        ax.hist(np.log10(spc + 1), bins=60, histtype="step", label=f"{eng} (n={len(spc)})")
    ax.set_xlabel("log10(sites covered per cell + 1)"); ax.set_ylabel("cells")
    ax.legend(); ax.set_title(f"{args.capture}: sites covered per cell")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "sites_per_cell.png"), dpi=120)
    plt.close(fig)

    # ---- pairwise matrix comparison: cellsnp vs souporcell (the headline) ----
    if "cellsnp" in engines and "souporcell" in engines:
        cs, cc, aA, tA, aB, tB = align(engines["cellsnp"], engines["souporcell"])
        report["cellsnp_vs_souporcell"] = {"common_sites": len(cs), "common_cells": len(cc)}

        # per-site total coverage
        sA = np.asarray(tA.sum(axis=1)).ravel(); sB = np.asarray(tB.sum(axis=1)).ravel()
        _loglog(sA, sB, "cellsnp site depth", "souporcell site depth",
                f"{args.capture}: per-site total depth", os.path.join(args.out, "site_depth.png"))
        report["cellsnp_vs_souporcell"]["site_depth_spearman"] = _spear(sA, sB)

        # per-cell total coverage
        cA = np.asarray(tA.sum(axis=0)).ravel(); cB = np.asarray(tB.sum(axis=0)).ravel()
        _loglog(cA, cB, "cellsnp cell depth", "souporcell cell depth",
                f"{args.capture}: per-cell total depth", os.path.join(args.out, "cell_depth.png"))
        report["cellsnp_vs_souporcell"]["cell_depth_spearman"] = _spear(cA, cB)

        # allele fraction at shared, well-covered (cell,site) entries
        tAc = tA.tocoo(); mask = {}
        # build dense AF comparison only over entries covered in both (memory-safe via COO intersect)
        A_alt = aA.tocsr(); A_tot = tA.tocsr(); B_alt = aB.tocsr(); B_tot = tB.tocsr()
        both = (A_tot >= args.min_cov).multiply(B_tot >= args.min_cov).tocoo()
        r, c = both.row, both.col
        if len(r):
            afA = np.asarray(A_alt[r, c]).ravel() / np.asarray(A_tot[r, c]).ravel()
            afB = np.asarray(B_alt[r, c]).ravel() / np.asarray(B_tot[r, c]).ravel()
            _hist2d(afA, afB, "cellsnp allele fraction", "souporcell allele fraction",
                    f"{args.capture}: AF concordance (cov>={args.min_cov}, n={len(r)})",
                    os.path.join(args.out, "af_concordance.png"))
            report["cellsnp_vs_souporcell"]["af_n_entries"] = int(len(r))
            report["cellsnp_vs_souporcell"]["af_median_abs_diff"] = float(np.median(np.abs(afA - afB)))
            report["cellsnp_vs_souporcell"]["af_frac_within_0.1"] = float(np.mean(np.abs(afA - afB) <= 0.1))
            report["cellsnp_vs_souporcell"]["af_pearson"] = float(np.corrcoef(afA, afB)[0, 1])

    with open(os.path.join(args.out, "metrics.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2))


def _spear(x, y):
    from scipy.stats import spearmanr
    m = (x > 0) | (y > 0)
    return float(spearmanr(x[m], y[m]).correlation) if m.sum() > 2 else float("nan")


def _loglog(x, y, xl, yl, title, path):
    m = (x > 0) & (y > 0)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.hexbin(np.log10(x[m]), np.log10(y[m]), gridsize=60, bins="log", cmap="viridis")
    lim = [min(np.log10(x[m]).min(), np.log10(y[m]).min()), max(np.log10(x[m]).max(), np.log10(y[m]).max())]
    ax.plot(lim, lim, "r--", lw=1)
    ax.set_xlabel(f"log10({xl})"); ax.set_ylabel(f"log10({yl})"); ax.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def _hist2d(x, y, xl, yl, title, path):
    fig, ax = plt.subplots(figsize=(5, 5))
    h = ax.hist2d(x, y, bins=50, range=[[0, 1], [0, 1]], cmap="magma", cmin=1)
    ax.plot([0, 1], [0, 1], "c--", lw=1); fig.colorbar(h[3], ax=ax, label="entries")
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


if __name__ == "__main__":
    main()
