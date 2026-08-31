import pandas as pd

# --- Sample sheet: one row per capture ------------------------------------
# Columns: capture_id, platform, bam, barcodes
captures = pd.read_csv(config["deps"]["captures"], dtype=str, index_col="capture_id")
captures.index = captures.index.astype(str)

_VALID_PLATFORMS = {"short", "long"}
_bad = set(captures["platform"].unique()) - _VALID_PLATFORMS
if _bad:
    # plain concatenation (see note in pileup_targets re: f-strings in .smk)
    raise ValueError("captures.csv has unknown platform(s) " + str(_bad)
                     + "; expected " + str(_VALID_PLATFORMS))

_VALID_ENGINES = {"cellsnp", "popscle", "souporcell"}
engines = config.get("engines", ["cellsnp"])
_bad_engines = set(engines) - _VALID_ENGINES
if _bad_engines:
    raise ValueError(f"config.engines has unknown engine(s) {_bad_engines}; "
                     f"valid engines: {sorted(_VALID_ENGINES)}")
_popscle_platforms = set(config["params"]["popscle"]["platforms"])
_souporcell_platforms = set(config["params"].get("souporcell", {}).get("platforms", []))


def lookup(field):
    """Return captures[field] for the wildcard capture_id."""
    def _inner(wildcards):
        return captures.loc[wildcards.capture_id, field]
    return _inner


def get_platform(wildcards):
    return captures.loc[wildcards.capture_id, "platform"]


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
        if "cellsnp" in engines:
            targets.append(str(out_dir / "cellsnp" / plat / capture_id / "cellSNP.tag.AD.mtx"))
        if "popscle" in engines and plat in _popscle_platforms:
            # request plp + umi explicitly: plp pulls cel/var (shared rules) but
            # not umi (produced by a separate rule).
            base = out_dir / "popscle" / plat / capture_id
            targets.append(str(base / "pileup.plp.gz"))
            targets.append(str(base / "pileup.umi.gz"))
        if "souporcell" in engines and plat in _souporcell_platforms:
            targets.append(str(out_dir / "souporcell" / plat / capture_id / "ref.mtx"))
            targets.append(str(out_dir / "souporcell" / plat / capture_id / "alt.mtx"))
    return targets


# --- Which VCF each capture piles at ----------------------------------------
# Three modes, resolved per capture by sites_vcf_input():
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
        raise ValueError("panel.build_in_workflow is set but captures.csv has no "
                         "'donors' column")
    val = captures.loc[wildcards.capture_id, "donors"]
    if pd.isna(val) or not str(val).strip():
        raise ValueError("panel.build_in_workflow is set but capture '"
                         + str(wildcards.capture_id) + "' has an empty 'donors' value")
    return str(val).strip()


def sites_vcf_input(wildcards):
    """The VCF the engines pile at for this capture (see the 3 modes above)."""
    if _panel_build:
        return str(out_dir / "panel" / wildcards.capture_id / "informative.vcf.gz")
    if "sites_vcf" in captures.columns:
        v = captures.loc[wildcards.capture_id, "sites_vcf"]
        if pd.notna(v) and str(v).strip():
            return str(v).strip()
    return config["deps"]["sites_vcf"]
