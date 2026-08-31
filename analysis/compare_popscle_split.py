#!/usr/bin/env python
"""Verify the chrom-split+merge popscle strategy against the unsplit ground truth.
Compares at the SEMANTIC level (barcode / genomic position), which is invariant
to the merge's droplet/SNP index relabeling.

Checks:
  1. var:  same SNP set (chrom:pos:ref:alt)
  2. cel:  per-barcode counts (NUM.READ/UMI/UMIwSNP/SNP) identical
  3. plp:  per-(barcode, SNP) allele pileup string identical
"""
import gzip
import sys

UNSPLIT = "/scratch/a56/jr9959/sc_pileup/results/popscle_unsplit/Neomet_Pool028a/unsplit"
SPLIT = "/scratch/a56/jr9959/sc_pileup/results/popscle/short/Neomet_Pool028a/pileup"


def op(p):
    return gzip.open(p, "rt")


def read_var(prefix):
    """SNP_ID(int) -> (chrom,pos,ref,alt) ; also the set of site keys."""
    d = {}
    with op(prefix + ".var.gz") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for l in fh:
            f = l.rstrip("\n").split("\t")
            sid = int(f[ix.get("#SNP_ID", ix.get("SNP_ID", 0))])
            key = (f[ix["CHROM"]], f[ix["POS"]], f[ix["REF"]], f[ix["ALT"]])
            d[sid] = key
    return d


def read_cel(prefix):
    """DROPLET_ID(int)->BARCODE ; BARCODE->counts tuple."""
    drop2bc, bc2cnt = {}, {}
    with op(prefix + ".cel.gz") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        did = ix.get("#DROPLET_ID", ix.get("DROPLET_ID"))
        cols = [c for c in ("NUM.READ", "NUM.UMI", "NUM.UMIwSNP", "NUM.SNP") if c in ix]
        for l in fh:
            f = l.rstrip("\n").split("\t")
            bc = f[ix["BARCODE"]]
            drop2bc[int(f[did])] = bc
            bc2cnt[bc] = tuple(f[ix[c]] for c in cols)
    return drop2bc, bc2cnt, cols


def read_plp(prefix, drop2bc, snp2key):
    """(barcode, site_key) -> ALLELES string (semantic key, index-invariant)."""
    d = {}
    with op(prefix + ".plp.gz") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        did = ix.get("#DROPLET_ID", ix.get("DROPLET_ID"))
        sid = ix.get("SNP_ID")
        al = ix.get("ALLELES")
        for l in fh:
            f = l.rstrip("\n").split("\t")
            bc = drop2bc.get(int(f[did])); key = snp2key.get(int(f[sid]))
            if bc is None or key is None:
                continue
            d[(bc, key)] = f[al] if al is not None else ""
    return d


print("reading var...")
uv, sv = read_var(UNSPLIT), read_var(SPLIT)
us, ss = set(uv.values()), set(sv.values())
print(f"1. VAR: unsplit={len(us)} split={len(ss)} shared={len(us & ss)} "
      f"only_unsplit={len(us - ss)} only_split={len(ss - us)}  "
      f"{'PASS' if us == ss else 'FAIL'}")

print("reading cel...")
ud, uc, cols = read_cel(UNSPLIT)
sd, sc, _ = read_cel(SPLIT)
bcs = set(uc) & set(sc)
mism = [b for b in bcs if uc[b] != sc[b]]
print(f"2. CEL ({','.join(cols)}): barcodes unsplit={len(uc)} split={len(sc)} shared={len(bcs)}  "
      f"count-mismatches={len(mism)}  {'PASS' if not mism and set(uc)==set(sc) else 'FAIL'}")
if mism[:3]:
    for b in mism[:3]:
        print(f"     e.g. {b}: unsplit={uc[b]} split={sc[b]}")

print("reading plp (this is the big one)...")
up = read_plp(UNSPLIT, ud, uv)
sp = read_plp(SPLIT, sd, sv)
keys = set(up) & set(sp)
pmis = sum(1 for k in keys if up[k] != sp[k])
print(f"3. PLP: entries unsplit={len(up)} split={len(sp)} shared={len(keys)}  "
      f"only_unsplit={len(set(up)-set(sp))} only_split={len(set(sp)-set(up))}  "
      f"allele-mismatches={pmis}  "
      f"{'PASS' if pmis==0 and set(up)==set(sp) else 'FAIL'}")

ok = (us == ss) and (not mism and set(uc) == set(sc)) and (pmis == 0 and set(up) == set(sp))
print(f"\n=== SPLIT-MERGE {'VERIFIED (identical to unsplit)' if ok else 'MISMATCH — investigate'} ===")
