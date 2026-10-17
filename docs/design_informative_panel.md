# Design — per-pool informative demux panels (tool + stage split)

## Goal

Genotype demux (vireo / demuxlet / freemuxlet / souporcell) needs a **donor
genotype VCF** restricted to the sites that (a) are reliably genotyped and (b)
actually distinguish the donors in *that* pool. Building that per-capture panel
is a small, pure transformation of an upstream cohort/project genotype VCF. This
doc fixes **where the logic lives, where the orchestration lives, and how the two
consumers (sc_pileup pileup, snp_demux assignment) stay in lock-step.**

The driving principle established during design:

> Separate three things that are easy to conflate — the **tool** (pure logic),
> the **workflow** (the pileup DAG, which must stay genotype-free and
> recyclable), and the **stage** (project-specific DVC orchestration).

## The three layers

```
  cohort / project germline VCF            (bcarc_wgs #57; or bcarc_snp imputed)
            │
            ▼  TOOL: tools/build_informative_panel.sh   (lives in sc_pileup, pure)
   per-capture panel/<capture>/informative.vcf.gz
            ▲
            │  STAGE: build_informative_panels           (lives in bcarc_chromium, DVC)
            │  supplies project pool-map + common-SNP backbone + genotype sources
            │
     ┌──────┴───────────────┐
     ▼                       ▼
  sc_pileup                snp_demux
  (sites_vcf)              (donor genotypes)
```

- **Tool** (`sc_pileup/tools/build_informative_panel.sh`) — reusable,
  project-agnostic. Ships with sc_pileup (the module that travels everywhere
  pileups are needed) but is **not a pileup rule**. Every project-specific input
  is an explicit argument; it reads no project config. One implementation, called
  from two entry points (below).
- **Workflow** (`sc_pileup/workflow`) — the pileup DAG stays genotype-free; it
  just piles at whatever `sites_vcf` it is handed. An optional thin rule
  (`rules/panel.smk`) calls the tool for standalone use, default off.
- **Stage** (bcarc_chromium `build_informative_panels`) — a DVC stage mirroring
  the existing `reheader_vcf` stage. It owns the project-specific facts (which
  cohort VCF, the captures→donors pool map from `config/snp_demux/`, the
  common-SNP backbone, the WGS/microarray source selection) and calls the tool in
  **batch** across all captures.

### Why the tool goes with sc_pileup but the stage goes with the project
- **Reuse across projects.** `#57` emits four project VCFs (brca_mega_atlas,
  neochemo_sc, phyllodes, luminal); each project's orchestrator calls the *same*
  tool — no copy-pasted bcftools logic per project. One implementation, no drift.
- **Single pool-map source.** The captures→donors map already lives
  authoritatively at the chromium level (`config/snp_demux/{captures,samples}.csv`,
  consumed by the snp_demux stage). The stage reads *that*; putting the split in
  sc_pileup would fork a second copy of pool composition.
- **sc_pileup stays genotype-free.** The one genotype-aware step doesn't belong
  in a generic pileup engine; as a sibling tool it doesn't pollute the DAG.
- **The stage is the natural home for source reconciliation.** Mixed WGS +
  microarray pools need per-donor source selection/merge — a project concern.

## The tool: composable steps, batch + single

`tools/build_informative_panel.sh` applies up to three optional steps, **folded
into one `bcftools view`** (so the informative allele-count is recomputed on the
donor subset):

| step | flag | effect |
|---|---|---|
| common-SNP intersect | `--common-sites REF` | restrict to a common-SNP backbone (`bcftools -T`); harmonises WGS vs microarray-imputed donors onto identical sites |
| pool subset | `--donors d1,d2,…` (required) | `bcftools -s` to the pool's donors |
| informative filter | `--informative` | `bcftools -c 1:minor`: keep only sites polymorphic among those donors (drop all-hom-ref / all-hom-alt) |

Modes: `--captures TSV` (batch; lines `capture_id⟶d1;d2;…`) or
`--capture ID --donors LIST` (single). Output: `<out-dir>/<capture>/informative.vcf.gz`.

**Purity rule:** the common-SNP reference, donor→source map, and pool map are all
*arguments* — never implicit config reads. This is what lets one tool serve every
project and every purpose.

### `--informative` vs `--common-sites` are different filters
They pull opposite ways on a **rare private variant** (present in one donor):
- `--informative` **keeps** it — it is maximally discriminative.
- `--common-sites` **drops** it — it is not population-common.

You usually want both for demux: the common-SNP backbone buys *reliability and
cross-platform harmonisation* (a rare variant seen only in the WGS donors is a
platform artifact, not signal), while the informative filter drops sites that
can't separate the pool. Known limitation: `-c 1:minor` still keeps the rare
all-heterozygous site (every donor 0/1); catching that needs a genotype-level
comparison and isn't a single bcftools flag — accepted as negligible for a
common-SNP panel.

### Purpose presets (set by the caller)
| purpose | `--common-sites` | `--informative` |
|---|---|---|
| demux | on | on |
| Track B (allelic imbalance, germline hets) | on | off |
| Track A (somatic) | off | off |

## The three pileup modes (sc_pileup side)

`sites_vcf_input()` (in `workflow/rules/common.smk`) picks per capture:

1. **`panel.build_in_workflow: true`** — sc_pileup builds the panel itself via
   `rule informative_panel` (needs a `donors` column in `captures.csv`). Standalone.
2. **per-capture `sites_vcf` column in `captures.csv`** — pile at a pre-built
   panel path. **Production demux path**: the chromium stage built the panels; we
   consume them.
3. **else `deps.sites_vcf`** — one global panel, or a somatic/Track-A VCF.

## Consumers stay in lock-step

The per-capture panel is **the same physical file** on both sides:
- **sc_pileup** consumes it as `sites_vcf` (what it piles at).
- **snp_demux** consumes it as the donor genotype VCF (`common.smk::panel_path()`
  → `vireo.smk::subset_vcf`, `split.smk::known_genotypes`).

One producer, two consumers ⇒ pileup sites ≡ demux genotypes by construction;
they can never silently drift and get intersected away.

## The common-SNP backbone

A **project-provided reference file** (pinned, DVC-tracked as a stage dep) — not
derived on the fly — so the backbone is reproducible and stamped into lineage.

Candidates investigated (GRCh38):

| candidate | ~#SNPs | threshold | notes |
|---|---|---|---|
| **souporcell `common_variants_grch38.vcf`** (IGSR 1000G native GRCh38) | ~10.4M | MAF≳0.02 | **recommended default**; native GRCh38 (no liftover); souporcell's built-in default so one file harmonises cellSNP/vireo/souporcell/popscle; distributed with [souporcell](https://github.com/wheaton5/souporcell) (pre-stage to a local path for offline/compute-node use) |
| cellSNP/vireo `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz` | ~7.4M | MAF≥0.05 | 2nd choice; the tool authors' documented list; hg19→hg38 liftover; 92 MB download ([SourceForge](https://sourceforge.net/projects/cellsnp/files/SNPlist/)) |
| gnomAD v3.1/v4 genomes, self-filtered | — | choose | best AFs, native GRCh38, but roll-your-own + heavy |
| dbSNP "common" | ~14M | MAF≥0.01 | **deprecated** — NCBI retired the common subset after build 151; avoid |
| Illumina GSA manifest | ~0.65M | array content | useful as a **QC overlay** to flag which backbone sites are directly-typed vs imputed for microarray donors; too sparse as the sole backbone |

**Recommended default:** the local native-GRCh38 souporcell 1000G common set.
Caveat before use: its records use **numeric chrom names (`1..X`, no `chr`)** —
match to the donor VCFs (`bcftools annotate --rename-chrs`) — and optionally
`bcftools view -m2 -M2 -v snps` to force strictly biallelic SNPs. Pin a
normalised copy as the project resource.

## Batch over per-capture (stage wiring)

The stage runs the tool in **batch** (one call, all captures) rather than a DVC
`foreach` over 100+ captures — a foreach bloats the lock file and status output.
Tradeoff: coarser DVC granularity (any capture/donor change reruns all panels),
but panel-building is cheap bcftools, so recompute cost is negligible and one
`outs: data/snp/vcf/panel/` directory is far simpler than 100 outs. A per-capture
mode remains in the tool for tests and the standalone rule.

## Migration status

Implemented in sc_pileup (this change):
- `tools/build_informative_panel.sh` (batch + single; the panel *rule* is now a
  thin wrapper around it, default off).
- `sites_vcf_input()` 3-mode routing; `panel.{build_in_workflow,common_sites,informative}`
  config; optional `sites_vcf` / `donors` columns in `captures.csv`.

Still to do (bcarc_chromium branch, upcoming):
- Add the `build_informative_panels` DVC stage (calls the tool in batch; deps =
  tool + project VCF + `config/snp_demux/{captures,samples}.csv` + common-SNP
  backbone; `outs: data/snp/vcf/panel/`).
- Repoint the snp_demux stage / sc_pileup stage to consume `data/snp/vcf/panel/`.
- First test: Atlas_Pool26a against `brca_mega_atlas.germline.common_snps.vcf.gz`
  + the pinned backbone.
