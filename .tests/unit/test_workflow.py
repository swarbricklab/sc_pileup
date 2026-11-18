#!/usr/bin/env python3
"""Small, data-free unit tests for sc_pileup rules.

Each test builds tiny synthetic fixtures (from the text files in data/) in a
scratch directory, writes a config, and runs Snakemake on specific targets with
the tools on PATH (no containers). Needs: snakemake 7.32, bcftools, samtools,
cellsnp-lite, pandas, pyyaml.

Usage (from the repo root):
    python .tests/unit/test_workflow.py            # all tests
    python .tests/unit/test_workflow.py somatic    # tests whose name matches
"""
import copy
import gzip
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
BASE_CONFIG = yaml.safe_load((HERE / "config.yaml").read_text())


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, text=True,
                          capture_output=True, **kw)


def deep_update(d, u):
    for k, v in u.items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            deep_update(d[k], v)
        else:
            d[k] = v
    return d


class Workdir:
    """A scratch dir holding fixtures, a captures.csv and a config."""

    def __init__(self, tmp: Path):
        self.root = tmp
        self.out = tmp / "out"
        self.logs = tmp / "logs"

    def path(self, name):
        return self.root / name

    def bgzip_vcf(self, src: Path, name: str) -> Path:
        dst = self.path(name)
        sh(f"bcftools view --no-version -Oz -o {dst} {src} && bcftools index -f {dst}")
        return dst

    def write_captures(self, rows, columns):
        p = self.path("captures.csv")
        lines = [",".join(columns)] + [",".join(str(r.get(c, "")) for c in columns) for r in rows]
        p.write_text("\n".join(lines) + "\n")
        return p

    def write_config(self, overrides=None, sites=None):
        cfg = copy.deepcopy(BASE_CONFIG)
        cfg["deps"]["captures"] = str(self.path("captures.csv"))
        cfg["deps"]["sites_vcf"] = str(sites or self.path("sites.vcf.gz"))
        cfg["outs"]["root"] = str(self.out)
        cfg["outs"]["logs"] = str(self.logs)
        cfg["ref"]["souporcell"]["genome_fasta"] = str(self.path("genome.fa"))
        deep_update(cfg, overrides or {})
        p = self.path("config.yaml")
        p.write_text(yaml.safe_dump(cfg, sort_keys=False))
        return p

    def profile(self):
        """The repo's workflow profile (set-resources etc.), minus its conda
        settings: the rules have no conda envs and the tools are on PATH, and
        use-conda makes snakemake call `conda info`, absent on CI runners."""
        d = self.path("profile")
        if not d.exists():
            prof = yaml.safe_load((REPO / "profiles/workflow/config.yaml").read_text())
            for k in ("use-conda", "conda-frontend"):
                prof.pop(k, None)
            d.mkdir()
            (d / "config.yaml").write_text(yaml.safe_dump(prof, sort_keys=False))
        return d

    def snakemake(self, targets=(), dry=False, check=True):
        cmd = [
            "snakemake", "--snakefile", str(REPO / "workflow/Snakefile"),
            "--configfile", str(self.path("config.yaml")),
            "--workflow-profile", str(self.profile()), "--retries", "0",
            "--directory", str(self.root), "--cores", "1",
            "--printshellcmds", "--nolock",
        ]
        if dry:
            cmd.append("-n")
        cmd += [str(t) for t in targets]
        res = subprocess.run(cmd, text=True, capture_output=True)
        if check and res.returncode != 0:
            raise AssertionError("snakemake failed:\n" + res.stdout + res.stderr)
        return res


def read_tsv(path):
    import pandas as pd
    return pd.read_csv(path, sep="\t", dtype=str)


def vcf_ids(path):
    with gzip.open(path, "rt") as fh:
        return [ln.split("\t")[2] for ln in fh if not ln.startswith("#")]


# --------------------------------------------------------------------------- #
# #16 somatic site panels
# --------------------------------------------------------------------------- #
def _somatic_workdir(w, captures_rows=None, columns=None):
    w.bgzip_vcf(DATA / "somatic_sites.vcf", "sites.vcf.gz")
    for f in ["pileup.bam", "pileup.bam.bai", "barcodes.tsv.gz", "genome.fa"]:
        w.path(f).touch()
    rows = captures_rows or [dict(capture_id="cap1", platform="short",
                                  bam=w.path("pileup.bam"),
                                  barcodes=w.path("barcodes.tsv.gz"))]
    w.write_captures(rows, columns or ["capture_id", "platform", "bam", "barcodes"])


def test_somatic_rejects_informative(w):
    _somatic_workdir(w)
    w.write_config({"panel_kind": "somatic", "panel": {"informative": True}})
    res = w.snakemake(dry=True, check=False)
    assert res.returncode != 0, "dry run should fail"
    out = res.stdout + res.stderr
    assert "germline-only panel filter(s) are set: panel.informative" in out, out


def test_somatic_rejects_common_sites_and_build(w):
    _somatic_workdir(w)
    w.write_config({"panel_kind": "somatic",
                    "panel": {"build_in_workflow": True, "common_sites": "x.vcf.gz"}})
    res = w.snakemake(dry=True, check=False)
    out = res.stdout + res.stderr
    assert res.returncode != 0
    assert "panel.build_in_workflow, panel.common_sites" in out, out


def test_somatic_warn_counts(w):
    _somatic_workdir(w)
    w.write_config({"panel_kind": "somatic"})
    report = w.out / "sites/cap1/site_types.tsv"
    w.snakemake([report])
    r = read_tsv(report).iloc[0].to_dict()
    assert r == {
        "capture_id": "cap1", "n_records": "5", "n_biallelic_snv": "3",
        "n_non_snv": "2", "n_multiallelic": "0", "n_mnv": "1", "n_indel": "1",
        "n_other": "0", "non_snv_action": "kept",
    }, r
    log = (w.logs / "sites/cap1/site_types.log").read_text()
    assert "WARNING: 2 of 5 record(s)" in log, log
    # warn mode: engines still pile at the original VCF, after the check
    dry = w.snakemake(dry=True).stdout
    assert "--regionsVCF " + str(w.path("sites.vcf.gz")) in dry, dry
    assert "biallelic_snv" not in dry


def test_somatic_drop(w):
    _somatic_workdir(w)
    w.write_config({"panel_kind": "somatic", "somatic": {"non_snv": "drop"}})
    sub = w.out / "sites/cap1/biallelic_snv.vcf.gz"
    w.snakemake([sub])
    assert vcf_ids(sub) == ["snv1", "snv2", "snv3"], vcf_ids(sub)
    r = read_tsv(w.out / "sites/cap1/site_types.tsv").iloc[0]
    assert r["non_snv_action"] == "dropped"
    dry = w.snakemake(dry=True).stdout
    assert "--regionsVCF " + str(sub) in dry, dry


def test_germline_capture_skips_site_check(w):
    cols = ["capture_id", "platform", "bam", "barcodes", "panel_kind"]
    rows = [
        dict(capture_id="germ", platform="short", bam=w.path("pileup.bam"),
             barcodes=w.path("barcodes.tsv.gz"), panel_kind=""),
        dict(capture_id="som", platform="short", bam=w.path("pileup.bam"),
             barcodes=w.path("barcodes.tsv.gz"), panel_kind="somatic"),
    ]
    _somatic_workdir(w, rows, cols)
    w.write_config()
    dry = w.snakemake(dry=True).stdout
    assert "sites/som/site_types.tsv" in dry
    assert "sites/germ/" not in dry


def test_bad_panel_kind(w):
    cols = ["capture_id", "platform", "bam", "barcodes", "panel_kind"]
    rows = [dict(capture_id="c", platform="short", bam=w.path("pileup.bam"),
                 barcodes=w.path("barcodes.tsv.gz"), panel_kind="tumour")]
    _somatic_workdir(w, rows, cols)
    w.write_config()
    res = w.snakemake(dry=True, check=False)
    assert res.returncode != 0 and "unknown panel_kind" in res.stdout + res.stderr


# --------------------------------------------------------------------------- #
# #17 nM/NH read pre-filter
# --------------------------------------------------------------------------- #
# data/reads.sam: 6 reads over chrT:100 (A>G), 2 barcodes, one UMI each.
#   r1 A nM0 NH1 | r2 G nM1 NH1 | r3 G nM7 NH1 | r4 A nM0 NH2 | r5 G nM2 NH1
#   r6 A, no nM/NH tags (kept: the filter only drops reads that fail it)
def _reads_workdir(w):
    sh(f"samtools view -b -o {w.path('pileup.bam')} {DATA / 'reads.sam'}"
       f" && samtools index {w.path('pileup.bam')}")
    sh(f"gzip -c {DATA / 'read_barcodes.tsv'} > {w.path('barcodes.tsv.gz')}")
    w.bgzip_vcf(DATA / "read_sites.vcf", "sites.vcf.gz")
    w.path("genome.fa").touch()
    w.write_captures([dict(capture_id="cap1", platform="short",
                           bam=w.path("pileup.bam"),
                           barcodes=w.path("barcodes.tsv.gz"))],
                     ["capture_id", "platform", "bam", "barcodes"])


def mtx_total(path):
    with open(path) as fh:
        body = [ln.split() for ln in fh if not ln.startswith("%")][1:]
    return sum(int(float(v)) for _, _, v in body)


def test_read_filter_off_is_noop(w):
    _reads_workdir(w)
    w.write_config({"params": {"read_filter": {"max_nM": None, "max_NH": None}}})
    dry = w.snakemake(dry=True).stdout
    assert "read_filter_bam" not in dry
    assert "--samFile " + str(w.path("pileup.bam")) in dry
    w.snakemake()
    cs = w.out / "cellsnp/short/cap1"
    assert mtx_total(cs / "cellSNP.tag.DP.mtx") == 6
    assert mtx_total(cs / "cellSNP.tag.AD.mtx") == 3


def test_read_filter_on_drops_nm_nh(w):
    _reads_workdir(w)
    w.write_config({"params": {"read_filter": {"max_nM": 5, "max_NH": 1}}})
    w.snakemake()
    cs = w.out / "cellsnp/short/cap1"
    # r3 (nM 7, ALT) and r4 (NH 2, REF) are gone: DP 6 -> 4, AD 3 -> 2
    assert mtx_total(cs / "cellSNP.tag.DP.mtx") == 4
    assert mtx_total(cs / "cellSNP.tag.AD.mtx") == 2
    log = (w.logs / "read_filter/short/cap1.log").read_text()
    assert "reads in:  6" in log and "reads out: 4" in log, log
    assert not (w.out / "read_filter/short/cap1/filtered.bam").exists()  # temp()


# --------------------------------------------------------------------------- #
def main(argv):
    pattern = argv[1] if len(argv) > 1 else ""
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f) and pattern in n]
    failed = []
    for name, fn in tests:
        tmp = Path(tempfile.mkdtemp(prefix="sc_pileup_unit_"))
        try:
            fn(Workdir(tmp))
            print("PASS", name)
            shutil.rmtree(tmp, ignore_errors=True)
        except Exception:
            failed.append(name)
            print("FAIL", name, "(scratch kept at", str(tmp) + ")")
            traceback.print_exc()
    print(f"\n{len(tests) - len(failed)}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
