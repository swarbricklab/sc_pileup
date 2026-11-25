#!/usr/bin/env python3
"""Small, data-free unit tests for sc_pileup rules.

Each test builds tiny synthetic fixtures (from the text files in data/) in a
scratch directory, writes a config, and runs Snakemake on specific targets with
the tools on PATH (no containers). Needs: snakemake 7.32, bcftools, samtools,
htslib (htsfile, bgzip), cellsnp-lite, pandas, pyyaml. bedtools is optional:
without it the old-pipe comparison of the covered-variant scan is skipped.

Usage (from the repo root):
    python .tests/unit/test_workflow.py            # all tests
    python .tests/unit/test_workflow.py somatic    # tests whose name matches
"""
import copy
import gzip
import os
import re
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

    def snakemake(self, targets=(), dry=False, check=True, extra=(), env=None):
        cmd = [
            "snakemake", "--snakefile", str(REPO / "workflow/Snakefile"),
            "--configfile", str(self.path("config.yaml")),
            "--workflow-profile", str(self.profile()), "--retries", "0",
            "--directory", str(self.root), "--cores", "1",
            "--printshellcmds", "--nolock",
        ]
        if dry:
            cmd.append("-n")
        cmd += list(extra)
        cmd += [str(t) for t in targets]
        res = subprocess.run(cmd, text=True, capture_output=True, env=env)
        if check and res.returncode != 0:
            raise AssertionError("snakemake failed:\n" + res.stdout + res.stderr)
        return res


class Skip(Exception):
    """Raised by a test that cannot run here (a tool is missing)."""


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


def test_somatic_all_snv_action_none(w):
    _somatic_workdir(w)
    snv_only = w.path("snv_only.vcf")
    snv_only.write_text("".join(
        l for l in (DATA / "somatic_sites.vcf").read_text().splitlines(True)
        if l.startswith("#") or "\tsnv" in l))
    w.bgzip_vcf(snv_only, "sites.vcf.gz")
    w.write_config({"panel_kind": "somatic", "somatic": {"non_snv": "drop"}})
    report = w.out / "sites/cap1/site_types.tsv"
    w.snakemake([report])
    r = read_tsv(report).iloc[0]
    assert (r["n_non_snv"], r["non_snv_action"]) == ("0", "none"), r.to_dict()


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
# Each read carries UR as well as UB: cellsnp-lite's default --UMItag Auto reads
# UR, and a fixture with UB only yields DP = 0 (real Cell Ranger BAMs have both).
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
    assert "reads lacking nM (kept): 1" in log, log
    assert "reads lacking NH (kept): 1" in log and "WARNING" not in log, log
    assert not (w.out / "read_filter/short/cap1/filtered.bam").exists()  # temp()


def test_read_filter_warns_when_tag_absent(w):
    _reads_workdir(w)
    # same reads with every nM tag stripped (as in e.g. a minimap2 BAM)
    sam = w.path("no_nm.sam")
    sam.write_text("".join(
        "\t".join(f for f in ln.rstrip("\n").split("\t")
                  if not f.startswith("nM:")) + "\n"
        for ln in (DATA / "reads.sam").read_text().splitlines(True)))
    sh(f"samtools view -b -o {w.path('pileup.bam')} {sam}"
       f" && samtools index {w.path('pileup.bam')}")
    w.write_config({"params": {"read_filter": {"max_nM": 5, "max_NH": None}}})
    res = w.snakemake()
    log = (w.logs / "read_filter/short/cap1.log").read_text()
    assert "reads lacking nM (kept): 6" in log and "NH" not in log, log
    assert "WARNING: no read in" in log and "nM cap removes nothing" in log, log
    assert "nM cap removes nothing" in res.stdout + res.stderr


# --------------------------------------------------------------------------- #
# #18 per-platform MAPQ / length
# --------------------------------------------------------------------------- #
def _cellsnp_cmds(dry_stdout):
    """{bam basename: 'cellsnp-lite ...' shell line} from a dry run with -p."""
    import re
    cmds = {}
    for m in re.finditer(r"cellsnp-lite --samFile (\S+).*", dry_stdout):
        cmds[Path(m.group(1)).name] = " ".join(m.group(0).split())
    return cmds


def _mixed_platform_workdir(w):
    for f in ["s.bam", "l.bam", "barcodes.tsv.gz", "sites.vcf.gz", "genome.fa"]:
        w.path(f).touch()
    w.write_captures(
        [dict(capture_id="cap_short", platform="short", bam=w.path("s.bam"),
              barcodes=w.path("barcodes.tsv.gz")),
         dict(capture_id="cap_long", platform="long", bam=w.path("l.bam"),
              barcodes=w.path("barcodes.tsv.gz"))],
        ["capture_id", "platform", "bam", "barcodes"])


def test_platform_params_mixed(w):
    _mixed_platform_workdir(w)
    w.write_config({"params": {"cellsnp_lite": {"platform_params": {
        "short": {}, "long": {"min_mapq": 7, "min_len": 200}}}}})
    cmds = _cellsnp_cmds(w.snakemake(dry=True).stdout)
    assert "--minMAPQ 20 --minLEN 30" in cmds["s.bam"], cmds
    assert "--minMAPQ 7 --minLEN 200" in cmds["l.bam"], cmds


def test_platform_params_fallback(w):
    _mixed_platform_workdir(w)
    # long overrides MAPQ only -> min_len falls back to the global 30;
    # no platform_params at all for short -> global values
    w.write_config({"params": {"cellsnp_lite": {"platform_params": {
        "long": {"min_mapq": 1}}}}})
    cmds = _cellsnp_cmds(w.snakemake(dry=True).stdout)
    assert "--minMAPQ 20 --minLEN 30" in cmds["s.bam"], cmds
    assert "--minMAPQ 1 --minLEN 30" in cmds["l.bam"], cmds


def test_platform_params_rejects_unknown_platform(w):
    _mixed_platform_workdir(w)
    w.write_config({"params": {"cellsnp_lite": {"platform_params": {
        "shrot": {"min_mapq": 1}}}}})
    res = w.snakemake(dry=True, check=False)
    out = res.stdout + res.stderr
    assert res.returncode != 0 and "shrot" in out, out


def test_shipped_config_long_placeholder(w):
    # the repo's config/config.yaml long-read placeholder is what long jobs get
    shipped = yaml.safe_load((REPO / "config/config.yaml").read_text())
    pp = shipped["params"]["cellsnp_lite"]["platform_params"]
    _mixed_platform_workdir(w)
    w.write_config({"params": {"cellsnp_lite": {"platform_params": pp}}})
    cmds = _cellsnp_cmds(w.snakemake(dry=True).stdout)
    assert "--minMAPQ 20 --minLEN 30" in cmds["s.bam"], cmds
    assert ("--minMAPQ %d --minLEN %d" % (pp["long"]["min_mapq"], pp["long"]["min_len"])
            in cmds["l.bam"]), cmds


# --------------------------------------------------------------------------- #
# #19 empty-droplet barcode selection
# --------------------------------------------------------------------------- #
# data/raw_matrix: 10 barcodes, 3 GEX features + 1 antibody feature.
# GEX UMI totals: AAAA 500 (cell), CCCC 50 (cell), GGGG 9, TTTT 10, ACAC 100,
# AGAG 101, ATAT 40 (+300 ADT), CACA 55, CGCG 0 (+20 ADT), CTCT 1.
RAW_GEX_TOTALS = {"AAAA-1": 500, "CCCC-1": 50, "GGGG-1": 9, "TTTT-1": 10,
                  "ACAC-1": 100, "AGAG-1": 101, "ATAT-1": 40, "CACA-1": 55,
                  "CGCG-1": 0, "CTCT-1": 1}
CELLS = ["AAAA-1", "CCCC-1"]


def _empties_workdir(w, raw_col=True):
    raw = w.path("raw_feature_bc_matrix")
    raw.mkdir()
    for f in ["matrix.mtx", "barcodes.tsv", "features.tsv"]:
        sh(f"gzip -c {DATA / 'raw_matrix' / f} > {raw / (f + '.gz')}")
    with gzip.open(w.path("barcodes.tsv.gz"), "wt") as fh:
        fh.write("".join(c + "\n" for c in CELLS))
    for f in ["pileup.bam", "pileup.bam.bai", "sites.vcf.gz", "genome.fa", "genome.fa.fai"]:
        w.path(f).touch()
    row = dict(capture_id="cap1", platform="short", bam=w.path("pileup.bam"),
               barcodes=w.path("barcodes.tsv.gz"), raw_matrix=raw)
    cols = ["capture_id", "platform", "bam", "barcodes"] + (["raw_matrix"] if raw_col else [])
    w.write_captures([row], cols)


def read_barcodes(path):
    with gzip.open(path, "rt") as fh:
        return [ln.rstrip("\n") for ln in fh]


def test_empties_selection(w):
    import pandas as pd
    _empties_workdir(w)
    w.write_config({"params": {"empties": {"umi_window": [10, 100]}}})
    bc = w.out / "empties/selection/cap1/barcodes.tsv.gz"
    w.snakemake([bc])
    got = read_barcodes(bc)
    assert got == ["ACAC-1", "ATAT-1", "CACA-1", "TTTT-1"], got
    assert not set(got) & set(CELLS)                       # disjoint from cells
    assert all(10 <= RAW_GEX_TOTALS[b] <= 100 for b in got)  # inside the window
    s = pd.read_csv(w.out / "empties/selection/cap1/summary.csv").iloc[0]
    assert s["n_selected"] == 4 and s["n_cell_barcodes_in_window"] == 1
    assert s["umi_window_min"] == 10 and s["umi_window_max"] == 100
    assert s["total_umis_selected"] == 205 and s["umi_median"] == 47.5, s
    # deterministic: a re-run gives byte-identical outputs
    first = bc.read_bytes(), (w.out / "empties/selection/cap1/summary.csv").read_bytes()
    w.snakemake([bc, "--forcerun", "select_empty_barcodes"])
    assert first == (bc.read_bytes(), (w.out / "empties/selection/cap1/summary.csv").read_bytes())


def test_empties_not_requested_by_default(w):
    _empties_workdir(w)
    w.write_config()
    dry = w.snakemake(dry=True).stdout
    assert "select_empty_barcodes" not in dry and "empties/" not in dry


def test_empties_needs_raw_matrix_column(w):
    _empties_workdir(w, raw_col=False)
    w.write_config()
    res = w.snakemake([w.out / "empties/selection/cap1/barcodes.tsv.gz"], dry=True, check=False)
    assert res.returncode != 0 and "no 'raw_matrix' column" in res.stdout + res.stderr


# --------------------------------------------------------------------------- #
# #20 barcode-set dimension (cells | empties)
# --------------------------------------------------------------------------- #
def _soup_bam(w, cbs):
    """A BAM with one 50 bp read over chrT:100 per (barcode, base) pair. Reads
    carry UR as well as UB: cellsnp-lite's --UMItag Auto reads UR."""
    seq = ("ACGT" * 12 + "AC")
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chrT\tLN:1000"]
    for i, (cb, base) in enumerate(cbs):
        sq = seq[:20] + base + seq[21:]
        lines.append("\t".join([f"r{i}", "0", "chrT", "80", "255", "50M", "*", "0",
                                 "0", sq, "I" * 50, f"CB:Z:{cb}",
                                 f"UB:Z:u{i}", f"UR:Z:u{i}"]))
    w.path("reads.sam").write_text("\n".join(lines) + "\n")
    sh(f"samtools view -b -o {w.path('pileup.bam')} {w.path('reads.sam')}"
       f" && samtools index {w.path('pileup.bam')}")


def test_barcode_sets_empties(w):
    _empties_workdir(w)
    w.path("pileup.bam").unlink()
    _soup_bam(w, [("AAAA-1", "A"), ("CCCC-1", "G"),       # cells
                  ("TTTT-1", "G"), ("CACA-1", "A"),       # empties (selected)
                  ("AGAG-1", "G")])                       # 101 UMIs: neither
    w.bgzip_vcf(DATA / "read_sites.vcf", "sites.vcf.gz")
    w.write_config({"barcode_sets": ["cells", "empties"]})
    w.snakemake()
    cells = w.out / "cellsnp/short/cap1"
    emp = w.out / "empties/cellsnp/short/cap1"
    # cells path unchanged; empties mirror it under empties/, nothing nested in cells
    assert (cells / "cellSNP.tag.AD.mtx").exists() and (emp / "cellSNP.tag.AD.mtx").exists()
    assert not any(p.is_dir() for p in cells.iterdir()), list(cells.iterdir())
    assert (w.logs / "empties/cellsnp/short/cap1.log").exists()
    selected = read_barcodes(w.out / "empties/selection/cap1/barcodes.tsv.gz")
    samples = (emp / "cellSNP.samples.tsv").read_text().split()
    assert samples == selected, (samples, selected)
    assert (cells / "cellSNP.samples.tsv").read_text().split() == CELLS
    assert mtx_total(emp / "cellSNP.tag.DP.mtx") == 2      # TTTT + CACA only
    assert mtx_total(emp / "cellSNP.tag.AD.mtx") == 1
    assert mtx_total(cells / "cellSNP.tag.DP.mtx") == 2


def test_barcode_sets_default_dag(w):
    _empties_workdir(w)
    w.write_config({"engines": ["cellsnp", "popscle", "souporcell"]})
    default = w.snakemake(dry=True).stdout
    assert "empties" not in default
    w.write_config({"engines": ["cellsnp", "popscle", "souporcell"],
                    "barcode_sets": ["cells", "empties"]})
    both = w.snakemake(dry=True).stdout
    assert "cellsnp_pileup_empties" in both and "select_empty_barcodes" in both
    # empties is cellsnp-only: no popscle / souporcell job touches it
    assert "empties/cellsnp/short/cap1/cellSNP.tag.AD.mtx" in both
    assert "empties/popscle" not in both and "empties/souporcell" not in both


def test_barcode_sets_validation(w):
    _empties_workdir(w)
    w.write_config({"barcode_sets": ["cells", "soup"]})
    res = w.snakemake(dry=True, check=False)
    assert res.returncode != 0
    w.write_config({"barcode_sets": ["empties"], "engines": ["popscle"]})
    res = w.snakemake(dry=True, check=False)
    assert res.returncode != 0 and "add 'cellsnp'" in res.stdout + res.stderr


# --------------------------------------------------------------------------- #
# #21 QC metrics for soup (empties) pileups
# --------------------------------------------------------------------------- #
OLD_METRIC_COLS = ["capture_id", "platform", "engine", "n_cells", "n_sites",
                   "n_sites_covered", "n_cells_covered", "total_depth", "total_alt",
                   "total_other", "alt_fraction", "mean_depth_per_covered_site"]


def test_qc_default_adds_only_barcode_set(w):
    import pandas as pd
    _reads_workdir(w)
    w.write_config({"qc": True})
    w.snakemake()
    for f in [w.out / "qc/short/cap1.metrics.csv", w.out / "qc/cohort_metrics.csv"]:
        m = pd.read_csv(f)
        expect = OLD_METRIC_COLS[:3] + ["barcode_set"] + OLD_METRIC_COLS[3:]
        assert list(m.columns) == expect, list(m.columns)
        assert m["barcode_set"].tolist() == ["cells"]
        r = m.iloc[0]
        assert (r["n_cells"], r["total_depth"], r["total_alt"]) == (2, 6, 3), r


def test_qc_empties_row(w):
    import pandas as pd
    _empties_workdir(w)
    w.path("pileup.bam").unlink()
    _soup_bam(w, [("AAAA-1", "A"), ("CCCC-1", "G"), ("TTTT-1", "G"),
                  ("CACA-1", "A"), ("AGAG-1", "G")])
    w.bgzip_vcf(DATA / "read_sites.vcf", "sites.vcf.gz")
    w.write_config({"qc": True, "barcode_sets": ["cells", "empties"]})
    w.snakemake()
    m = pd.read_csv(w.out / "qc/cohort_metrics.csv")
    assert m["barcode_set"].tolist() == ["cells", "empties"], m
    assert m["engine"].tolist() == ["cellsnp", "cellsnp"]
    cells, emp = m.iloc[0], m.iloc[1]
    assert pd.isna(cells["umis_per_barcode_median"])
    assert emp["n_cells"] == 4                 # selected empty-droplet barcodes
    assert emp["n_cells_covered"] == 2 and emp["total_depth"] == 2
    assert emp["total_alt"] == 1 and emp["alt_fraction"] == 0.5
    assert emp["umis_per_barcode_median"] == 47.5


# --------------------------------------------------------------------------- #
# souporcell covered-variant scan (scripts/covered_variants.sh)
# --------------------------------------------------------------------------- #
# data/scan_reads.sam gives chrT these read depths (0 everywhere else):
#   100-104: 2 | 105: 0 | 106-108: 3 | 109: 4 | 110: 1 | 111-112: 2
# With min_cov 2 and max_cov 4, the in-range runs are 100-104, 106-108 and
# 111-112. A record at POS with a REF of length L is tested at depth positions
# POS-1 .. POS+L-2 and written once for each in-range run among them. Some
# records of data/scan_sites.vcf:
#   s_own_in (POS 100) is tested at 99, depth 0: dropped
#   s_before_in (101) is tested at 100, depth 2: kept
#   del_two_runs (104, 6 bases) is tested at 103-108, two runs: written twice
#   del_three_runs (104, 10 bases) is tested at 103-112, three runs: three times
#   s_at_min (105) is tested at 104, depth 2: kept
#   s_at_max (110) is tested at 109, depth 4: dropped
#   s_below_min (111) is tested at 110, depth 1: dropped
#   pos_one (1) is tested at no position, and chrU has no reads: both dropped
SCAN_KEPT = ["s_before_in", "del_one_run", "s_run_middle",
             "del_two_runs", "del_two_runs",
             "del_three_runs", "del_three_runs", "del_three_runs",
             "s_at_min", "gap_into_run", "ins", "mnv", "multi"]
SCAN_OUT = "souporcell/short/cap1/covered_variants.vcf"

# The rule's shell before scripts/covered_variants.sh (souporcell's whole-BAM
# pipe), with @BAM@, @VCF@, @LO@, @HI@, @BED@ and @OUT@ for its inputs.
OLD_SCAN = r"""
samtools depth -d 0 @BAM@ \
    | awk -v lo=@LO@ -v hi=@HI@ '{ if ($3 >= lo && $3 < hi) print $1"\t"$2"\t"($2+1)"\t"$3 }' \
    | bedtools merge -i - >@BED@
(
    bcftools view -h @VCF@ | grep -v '^##bcftools_'
    bedtools intersect -wa -a @VCF@ -b @BED@
) >@OUT@
"""


def _scan_workdir(w):
    sh(f"samtools view -b -o {w.path('pileup.bam')} {DATA / 'scan_reads.sam'}"
       f" && samtools index {w.path('pileup.bam')}")
    for f in ["barcodes.tsv.gz", "genome.fa"]:
        w.path(f).touch()
    w.write_captures([dict(capture_id="cap1", platform="short",
                           bam=w.path("pileup.bam"),
                           barcodes=w.path("barcodes.tsv.gz"))],
                     ["capture_id", "platform", "bam", "barcodes"])


def _scan(w, sites, check=True):
    """Run souporcell_covered_variants on `sites` (min_cov 2, max_cov 4)."""
    w.write_config({"engines": ["souporcell"],
                    "params": {"souporcell": {"min_cov": 2, "max_cov": 4}}},
                   sites=sites)
    out = w.out / SCAN_OUT
    if out.exists():
        out.unlink()
    res = w.snakemake([out], check=check)
    return out, res


def _vcf_parts(path):
    lines = Path(path).read_text().splitlines()
    header = [ln for ln in lines if ln.startswith("#")]
    ids = [ln.split("\t")[2] for ln in lines if not ln.startswith("#")]
    return header, ids


def test_covered_variants_records(w):
    _scan_workdir(w)
    sites = w.bgzip_vcf(DATA / "scan_sites.vcf", "sites.vcf.gz")
    out, _ = _scan(w, sites)
    header, ids = _vcf_parts(out)
    assert ids == SCAN_KEPT, ids
    assert not [h for h in header if h.startswith("##bcftools_")], header


def test_covered_variants_plain_vcf(w):
    _scan_workdir(w)
    sites = w.path("sites.vcf")
    shutil.copy(DATA / "scan_sites.vcf", sites)
    out, _ = _scan(w, sites)
    header, ids = _vcf_parts(out)
    assert ids == SCAN_KEPT, ids
    # htslib's header: PASS line added, the ##bcftools_ line dropped
    want = sh(f"bcftools view -h {sites} | grep -v '^##bcftools_'").stdout
    assert header == want.splitlines(), (header, want)
    assert '##FILTER=<ID=PASS,Description="All filters passed">' in header


def test_covered_variants_matches_old_pipe(w):
    if not shutil.which("bedtools"):
        raise Skip("bedtools is not on PATH")
    _scan_workdir(w)
    gz = w.bgzip_vcf(DATA / "scan_sites.vcf", "sites.vcf.gz")
    plain = w.path("sites.vcf")
    shutil.copy(DATA / "scan_sites.vcf", plain)
    for sites in (gz, plain):
        want = w.path(sites.name + ".old_pipe.vcf")
        cmd = OLD_SCAN
        for k, v in [("@BAM@", w.path("pileup.bam")), ("@VCF@", sites),
                     ("@LO@", 2), ("@HI@", 4), ("@OUT@", want),
                     ("@BED@", w.path("depth_merged.bed"))]:
            cmd = cmd.replace(k, str(v))
        subprocess.run(["bash", "-c", "set -euo pipefail\n" + cmd], check=True,
                       text=True, capture_output=True)
        out, _ = _scan(w, sites)
        assert out.read_bytes() == want.read_bytes(), (
            sites.name, out.read_text(), want.read_text())


def test_scan_container_default_and_override(w):
    # The scan runs in containers.samtools_depth, or in samtools 1.22 when the
    # config has no such key; split_bam_by_chromosome keeps containers.samtools.
    # A dry run with --use-singularity names the image each requested job would
    # pull. A stand-in `singularity` that only reports a version lets that dry
    # run work where singularity is not installed.
    _scan_workdir(w)
    w.path("sites.vcf.gz").touch()
    shim = w.path("bin")
    shim.mkdir()
    (shim / "singularity").write_text("#!/bin/sh\necho 'singularity version 3.8.0'\n")
    (shim / "singularity").chmod(0o755)
    env = dict(os.environ, PATH=str(shim) + os.pathsep + os.environ["PATH"])
    sing = ["--use-singularity", "--singularity-prefix", str(w.path("images"))]

    def images(target):
        res = w.snakemake([target], dry=True, extra=sing, env=env)
        return re.findall(r"Singularity image (\S+) will be pulled", res.stdout + res.stderr)

    scan = w.out / SCAN_OUT
    split = w.out / "popscle/short/cap1/split/chr1.bam"
    default = "docker://quay.io/biocontainers/samtools:1.22--h96c455f_0"
    shared = BASE_CONFIG["containers"]["samtools"]
    assert "samtools_depth" not in BASE_CONFIG["containers"]
    w.write_config({"engines": ["souporcell"]})
    assert images(scan) == [default], images(scan)
    assert images(split) == [shared], images(split)
    other = "docker://example.org/samtools:test"
    w.write_config({"engines": ["souporcell"],
                    "containers": {"samtools_depth": other}})
    assert images(scan) == [other], images(scan)
    assert images(split) == [shared], images(split)


def test_covered_variants_symbolic_alt_fails(w):
    _scan_workdir(w)
    sites = w.path("sites.vcf")
    shutil.copy(DATA / "scan_symbolic.vcf", sites)
    _, res = _scan(w, sites, check=False)
    assert res.returncode != 0, "the scan should refuse a symbolic ALT"
    log = (w.logs / "souporcell/short/cap1/covered_variants.log").read_text()
    assert "symbolic ALT or POS < 1" in log, log


# --------------------------------------------------------------------------- #
def main(argv):
    pattern = argv[1] if len(argv) > 1 else ""
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f) and pattern in n]
    failed, skipped = [], []
    for name, fn in tests:
        tmp = Path(tempfile.mkdtemp(prefix="sc_pileup_unit_"))
        try:
            fn(Workdir(tmp))
            print("PASS", name)
            shutil.rmtree(tmp, ignore_errors=True)
        except Skip as why:
            skipped.append(name)
            print("SKIP", name, f"({why})")
            shutil.rmtree(tmp, ignore_errors=True)
        except Exception:
            failed.append(name)
            print("FAIL", name, "(scratch kept at", str(tmp) + ")")
            traceback.print_exc()
    passed = len(tests) - len(failed) - len(skipped)
    print(f"\n{passed}/{len(tests)} passed, {len(skipped)} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
