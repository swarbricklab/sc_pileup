import pandas as pd

# --- Sample sheet: one row per capture ------------------------------------
# Columns: capture_id, platform, bam, barcodes
captures = pd.read_csv(config["deps"]["captures"], dtype=str, index_col="capture_id")
captures.index = captures.index.astype(str)

_VALID_PLATFORMS = {"short", "long"}
_bad = set(captures["platform"].unique()) - _VALID_PLATFORMS
if _bad:
    # plain concatenation (see note in pileup_targets re: f-strings in .smk)
    raise ValueError(
        "captures.csv has unknown platform(s) "
        + str(_bad)
        + "; expected "
        + str(_VALID_PLATFORMS)
    )

_VALID_ENGINES = {"cellsnp", "popscle", "souporcell"}
engines = config.get("engines", ["cellsnp"])
_bad_engines = set(engines) - _VALID_ENGINES
if _bad_engines:
    raise ValueError(
        f"config.engines has unknown engine(s) {_bad_engines}; "
        f"valid engines: {sorted(_VALID_ENGINES)}"
    )
# Barcode sets each capture is piled up over. `cells` (the cell-called
# `barcodes` column) keeps today's paths; `empties` (empty droplets from
# rules/empties.smk, for soup genotyping) runs cellsnp only and lands one level
# down, in results/cellsnp/<platform>/<capture_id>/empties/.
_VALID_BARCODE_SETS = {"cells", "empties"}
barcode_sets = config.get("barcode_sets", ["cells"])
_bad_sets = set(barcode_sets) - _VALID_BARCODE_SETS
if _bad_sets:
    raise ValueError(
        "config.barcode_sets has unknown set(s) "
        + str(sorted(_bad_sets))
        + "; valid: "
        + str(sorted(_VALID_BARCODE_SETS))
    )
if "empties" in barcode_sets:
    if "cellsnp" not in engines:
        raise ValueError(
            "barcode_sets includes 'empties', which is piled up by the cellsnp "
            "engine only; add 'cellsnp' to config.engines"
        )
    if "raw_matrix" not in captures.columns:
        raise ValueError(
            "barcode_sets includes 'empties' but captures.csv has no 'raw_matrix' "
            "column (path to each capture's Cell Ranger raw_feature_bc_matrix/)"
        )


# capture_id is one path component; this keeps .../<capture_id>/empties/... from
# also matching the cells rules with capture_id="<id>/empties".
wildcard_constraints:
    capture_id="[^/]+",


_popscle_platforms = set(config["params"]["popscle"]["platforms"])
_souporcell_platforms = set(config["params"].get("souporcell", {}).get("platforms", []))


def lookup(field):
    """Return captures[field] for the wildcard capture_id."""

    def _inner(wildcards):
        return captures.loc[wildcards.capture_id, field]

    return _inner


def get_platform(wildcards):
    return captures.loc[wildcards.capture_id, "platform"]


# --- Optional nM/NH read pre-filter (cellsnp engine; see rules/read_filter.smk)
_read_filter = config["params"].get("read_filter") or {}
_read_filter_terms = []
read_filter_tags = []  # the tags actually capped (logged: how many reads lack each)
for _tag, _key in [("nM", "max_nM"), ("NH", "max_NH")]:
    if _read_filter.get(_key) is not None:
        read_filter_tags.append(_tag)
        # keep reads lacking the tag; drop only those that exceed the cap
        _read_filter_terms.append(
            "(!["
            + _tag
            + "] || ["
            + _tag
            + "] <= "
            + str(int(_read_filter[_key]))
            + ")"
        )
read_filter_expr = " && ".join(_read_filter_terms)


def cellsnp_bam(wildcards):
    """The BAM cellsnp reads: the nM/NH-filtered copy when params.read_filter is
    on, else the capture's own BAM."""
    if read_filter_expr:
        return str(
            out_dir
            / "read_filter"
            / get_platform(wildcards)
            / wildcards.capture_id
            / "filtered.bam"
        )
    return captures.loc[wildcards.capture_id, "bam"]


def cellsnp_bai(wildcards):
    """The filtered BAM's index (so temp() keeps it until cellsnp is done);
    [] when the filter is off, leaving the original cellsnp inputs unchanged."""
    return (cellsnp_bam(wildcards) + ".bai") if read_filter_expr else []


def cellsnp_platform_param(name):
    """Input-function for a per-platform cellsnp read-admission threshold:
    params.cellsnp_lite.platform_params.<platform>.<name>, falling back to the
    global params.cellsnp_lite.<name>."""
    cs = config["params"]["cellsnp_lite"]

    def _inner(wildcards):
        plat = captures.loc[wildcards.capture_id, "platform"]
        per = (cs.get("platform_params") or {}).get(plat) or {}
        v = per.get(name)
        return cs[name] if v is None else v

    return _inner


def cellsnp_extra(wildcards):
    """Per-platform extra cellsnp-lite flags."""
    plat = captures.loc[wildcards.capture_id, "platform"]
    return config["params"]["cellsnp_lite"].get("platform_args", {}).get(plat, "")


def pileup_targets():
    """All requested pileup outputs across captures × engines (platform-aware)."""
    # NB: build paths with pathlib `/`-chaining, NOT f-strings. Snakemake's
    # source reconstruction inserts spaces around `{...}` in f-strings that live
    # in .smk function bodies, which corrupts the resulting paths.
    targets = []
    for capture_id, row in captures.iterrows():
        plat = row["platform"]
        if "empties" in barcode_sets:
            targets.append(
                str(
                    out_dir
                    / "cellsnp"
                    / plat
                    / capture_id
                    / "empties"
                    / "cellSNP.tag.AD.mtx"
                )
            )
        if "cells" not in barcode_sets:
            continue
        if "cellsnp" in engines:
            targets.append(
                str(out_dir / "cellsnp" / plat / capture_id / "cellSNP.tag.AD.mtx")
            )
        if "popscle" in engines and plat in _popscle_platforms:
            # request plp + umi explicitly: plp pulls cel/var (shared rules) but
            # not umi (produced by a separate rule).
            base = out_dir / "popscle" / plat / capture_id
            targets.append(str(base / "pileup.plp.gz"))
            targets.append(str(base / "pileup.umi.gz"))
        if "souporcell" in engines and plat in _souporcell_platforms:
            targets.append(str(out_dir / "souporcell" / plat / capture_id / "ref.mtx"))
            targets.append(str(out_dir / "souporcell" / plat / capture_id / "alt.mtx"))
    # Additive QC summary (goal-1: reads the outputs above, changes nothing).
    if config.get("qc", True):
        targets.append(str(out_dir / "qc/cohort_metrics.csv"))
    return targets


# --- Which VCF each capture piles at ----------------------------------------
# Three modes, resolved per capture by raw_sites_vcf() (sites_vcf_input() adds
# the optional somatic biallelic-SNV subset on top -- see "Panel kind" below):
#   1. panel.build_in_workflow: sc_pileup builds the per-pool panel itself via
#      rule informative_panel (tools/build_informative_panel.sh) -> pile at it.
#      Standalone convenience; needs a `donors` column in captures.csv.
#   2. per-capture `sites_vcf` column in captures.csv: pile at that path. This is
#      the production path -- an upstream orchestration stage (e.g. bcarc_chromium
#      build_informative_panels) built the panels; we just consume them.
#   3. else: the global deps.sites_vcf (one shared panel, or a somatic/Track-A VCF).
_panel = config.get("panel", {})
_panel_build = bool(_panel.get("build_in_workflow", False))
_panel_common = (_panel.get("common_sites") or "").strip()
_panel_informative = bool(_panel.get("informative", True))
# Arg strings for the tool (module-level: avoid f-strings inside defs -- see note).
panel_common_arg = ("--common-sites " + _panel_common) if _panel_common else ""
panel_informative_arg = "--informative" if _panel_informative else ""


def pool_donors(wildcards):
    """';'/','-separated donor sample names for this capture (for the tool's
    --donors / bcftools -s). Only needed when panel.build_in_workflow is set."""
    if "donors" not in captures.columns:
        raise ValueError(
            "panel.build_in_workflow is set but captures.csv has no " "'donors' column"
        )
    val = captures.loc[wildcards.capture_id, "donors"]
    if pd.isna(val) or not str(val).strip():
        raise ValueError(
            "panel.build_in_workflow is set but capture '"
            + str(wildcards.capture_id)
            + "' has an empty 'donors' value"
        )
    return str(val).strip()


def raw_sites_vcf(wildcards):
    """The site VCF resolved for this capture (see the 3 modes above), before
    any somatic site-type filtering."""
    if _panel_build:
        return str(out_dir / "panel" / wildcards.capture_id / "informative.vcf.gz")
    if "sites_vcf" in captures.columns:
        v = captures.loc[wildcards.capture_id, "sites_vcf"]
        if pd.notna(v) and str(v).strip():
            return str(v).strip()
    return config["deps"]["sites_vcf"]


def sites_vcf_input(wildcards):
    """The VCF the engines pile at for this capture: the resolved site VCF, or
    (somatic capture with somatic.non_snv: drop) its biallelic-SNV subset."""
    if panel_kind(wildcards.capture_id) == "somatic" and _somatic_non_snv == "drop":
        return str(out_dir / "sites" / wildcards.capture_id / "biallelic_snv.vcf.gz")
    return raw_sites_vcf(wildcards)


# --- Panel kind: germline (default) | somatic --------------------------------
# A somatic (Track A) site VCF must not go through the germline panel filters:
# `panel.informative` keeps only sites polymorphic among the pool's genotypes and
# `panel.common_sites` restricts to a common-SNP backbone -- either would silently
# delete most somatic sites. So for somatic captures we refuse those settings at
# DAG build, and check the site VCF's record types before piling up (cellSNP-lite
# piles up biallelic SNVs; indels/MNVs only inflate the "unreachable" fraction).
# Set globally via `panel_kind`, overridable per capture with a `panel_kind`
# column in captures.csv (blank = global value).
_VALID_PANEL_KINDS = {"germline", "somatic"}
_VALID_NON_SNV = {"warn", "drop"}
_global_panel_kind = config.get("panel_kind", "germline")
_somatic_non_snv = config.get("somatic", {}).get("non_snv", "warn")
if _somatic_non_snv not in _VALID_NON_SNV:
    raise ValueError(
        "somatic.non_snv must be one of "
        + str(sorted(_VALID_NON_SNV))
        + ", got '"
        + str(_somatic_non_snv)
        + "'"
    )


def panel_kind(capture_id):
    """germline | somatic for this capture (captures.csv column, else global)."""
    if "panel_kind" in captures.columns:
        v = captures.loc[capture_id, "panel_kind"]
        if pd.notna(v) and str(v).strip():
            return str(v).strip()
    return _global_panel_kind


_bad_kinds = {
    str(c) + "=" + panel_kind(c)
    for c in captures.index
    if panel_kind(c) not in _VALID_PANEL_KINDS
}
if _bad_kinds:
    raise ValueError(
        "unknown panel_kind "
        + str(sorted(_bad_kinds))
        + "; expected one of "
        + str(sorted(_VALID_PANEL_KINDS))
    )
somatic_captures = [c for c in captures.index if panel_kind(c) == "somatic"]
_germline_only = [
    name
    for name, on in [
        ("panel.build_in_workflow", _panel_build),
        ("panel.informative", _panel_informative),
        ("panel.common_sites", bool(_panel_common)),
    ]
    if on
]
if somatic_captures and _germline_only:
    raise ValueError(
        "panel_kind is 'somatic' for capture(s) "
        + str(somatic_captures)
        + " but germline-only panel filter(s) are set: "
        + ", ".join(_germline_only)
        + ". These would silently delete most somatic sites. For a somatic "
        "(Track A) site VCF set panel.build_in_workflow: false, "
        "panel.informative: false and panel.common_sites: null."
    )


def site_check_input(wildcards):
    """The somatic site-type report, for somatic captures only ([] otherwise).
    Engines take it as an input so the record-type check (and its warning) runs
    before any pileup."""
    if panel_kind(wildcards.capture_id) == "somatic":
        return str(out_dir / "sites" / wildcards.capture_id / "site_types.tsv")
    return []
