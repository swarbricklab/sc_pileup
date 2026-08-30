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

engines = config.get("engines", ["cellsnp"])
_popscle_platforms = set(config["params"]["popscle"]["platforms"])


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
    return targets
