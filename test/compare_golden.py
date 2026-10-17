#!/usr/bin/env python3
"""Compare a fresh sc_pileup test run against the DVC-tracked golden.

The pileup payload is deterministic, but gzip wrappers embed an mtime and VCF
headers carry volatile `##` lines, so we compare *content*, not raw bytes:

  * .vcf / .vcf.gz : compare the data records only (lines not starting with '#')
  * other .gz      : compare the DECOMPRESSED bytes (ignores the gzip mtime)
  * everything else: compare raw bytes

Exits non-zero (and lists them) if any golden file is missing from the fresh
output or differs in content. Index files are ignored.

Usage:  compare_golden.py GOLDEN_DIR FRESH_DIR
"""
import gzip
import hashlib
import sys
from pathlib import Path

IGNORE_SUFFIXES = (".bai", ".csi", ".tbi", ".fai")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def content_key(path: Path) -> str:
    """A content fingerprint that ignores gzip mtime and VCF header volatility."""
    name = path.name
    if name.endswith(".vcf") or name.endswith(".vcf.gz"):
        opener = gzip.open if name.endswith(".gz") else open
        with opener(path, "rt") as fh:
            records = [ln for ln in fh if not ln.startswith("#")]
        return _digest("".join(records).encode())
    if name.endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return _digest(fh.read())
    return _digest(path.read_bytes())


def collect(root: Path):
    return {
        p.relative_to(root).as_posix(): p
        for p in sorted(root.rglob("*"))
        if p.is_file() and not p.name.endswith(IGNORE_SUFFIXES)
    }


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    golden_dir, fresh_dir = Path(sys.argv[1]), Path(sys.argv[2])
    golden, fresh = collect(golden_dir), collect(fresh_dir)

    missing, differ = [], []
    for rel, gpath in golden.items():
        fpath = fresh.get(rel)
        if fpath is None:
            missing.append(rel)
        elif content_key(gpath) != content_key(fpath):
            differ.append(rel)

    extra = sorted(set(fresh) - set(golden))

    print(f"golden files: {len(golden)}   fresh files: {len(fresh)}")
    print(f"OK:      {len(golden) - len(missing) - len(differ)}")
    for rel in missing:
        print(f"MISSING: {rel}")
    for rel in differ:
        print(f"DIFFER:  {rel}")
    for rel in extra:
        print(f"extra (ignored): {rel}")

    if missing or differ:
        print("\nFAIL: fresh output does not match the golden.")
        return 1
    print("\nPASS: fresh output matches the golden.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
