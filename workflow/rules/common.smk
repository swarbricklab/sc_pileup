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


# --- Optional per-pool informative panel ------------------------------------
# When config.panel.informative_filter is true, engines pile at a per-capture
# panel (sites_vcf subset to the pool's donors, kept only where polymorphic among
# them; see rules/panel.smk) instead of the raw sites_vcf. Correct for demux;
# leave FALSE for somatic / Track-A pileups (no donor subsetting wanted).
_informative = bool(config.get("panel", {}).get("informative_filter", False))


def pool_donors(wildcards):
    """';'-separated donor sample names for this capture (for bcftools -s)."""
    if "donors" not in captures.columns:
        raise ValueError("config.panel.informative_filter is set but captures.csv "
                         "has no 'donors' column")
    val = captures.loc[wildcards.capture_id, "donors"]
    if pd.isna(val) or not str(val).strip():
        raise ValueError("config.panel.informative_filter is set but capture '"
                         + str(wildcards.capture_id) + "' has an empty 'donors' value")
    return str(val).strip()


def sites_vcf_input(wildcards):
    """VCF the engines pile at: the per-pool informative panel when
    informative_filter is on, else the global sites_vcf. Built with pathlib to
    match rule informative_panel's output string (see f-string note above)."""
    if _informative:
        return str(out_dir / "panel" / wildcards.capture_id / "informative.vcf.gz")
    return config["deps"]["sites_vcf"]
