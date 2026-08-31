#!/usr/bin/env bash
# merge_genotype_sources.sh -- combine per-source donor genotypes into ONE
# multi-sample VCF for a MIXED-source pool (e.g. WGS + microarray-imputed).
#
# Reusable, project-agnostic (ships with sc_pileup; called by orchestration
# stages). Each --source is "<vcf>:<donor1,donor2,...>": the tool subsets that
# VCF to those donors, then bcftools-merges all sources at the union of sites
# (samples = all donors; a donor is missing where its source lacks a site),
# normalises, and optionally fills missing GTs to hom-ref. Sample names in each
# source must already be donor IDs (so the merged columns are donor IDs).
#
#   --out FILE            output bgzipped VCF (indexed)
#   --source VCF:DONORS   repeatable; per-source donor subset (','-separated)
#   [--regions FILE]      restrict every source to these sites first (BED/tab/VCF)
#   [--fill-ref]          set missing GTs to 0/0 after merge (bcftools +setGT).
#                         Reasonable for WGS (no call = hom-ref); for microarray
#                         donors it assumes unmeasured = hom-ref, so keep the
#                         merge restricted to sites both sources genotype.
#
# The mechanism is generic; the donor->source POLICY (which donor from which
# VCF) is the caller's, passed as --source args. Requires bcftools (+ setGT
# plugin if --fill-ref); run in the bcftools/demuxafy container.
set -euo pipefail

out=""; regions=""; fill=0; srcs=()
while [ $# -gt 0 ]; do case "$1" in
  --out)     out=$2; shift 2;;
  --source)  srcs+=("$2"); shift 2;;
  --regions) regions=$2; shift 2;;
  --fill-ref) fill=1; shift;;
  -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
  *) echo "unknown arg: $1" >&2; exit 1;;
esac; done
[ -n "$out" ] && [ ${#srcs[@]} -ge 1 ] || { echo "need --out and >=1 --source" >&2; exit 1; }

tmpd=$(mktemp -d); trap 'rm -rf "$tmpd"' EXIT
parts=(); i=0
for s in "${srcs[@]}"; do
  vcf=${s%%:*}; donors=${s#*:}; donors=${donors//;/,}
  p="$tmpd/src$i.vcf.gz"
  args=(view -s "$donors" --force-samples)
  [ -n "$regions" ] && args+=(-T "$regions")
  bcftools "${args[@]}" "$vcf" -Oz -o "$p"
  bcftools index "$p"
  parts+=("$p"); i=$((i+1))
done

if [ ${#parts[@]} -eq 1 ]; then
  merged="${parts[0]}"
else
  merged="$tmpd/merged.vcf.gz"
  bcftools merge "${parts[@]}" -Oz -o "$merged"
  bcftools index "$merged"
fi

if [ "$fill" = 1 ]; then
  # missing GT -> hom-ref. setGT is a standard plugin; fail loudly if absent.
  bcftools +setGT "$merged" -Oz -o "$out" -- -t . -n 0
else
  bcftools view "$merged" -Oz -o "$out"
fi
bcftools index -f "$out"
echo "merged $(bcftools query -l "$out" | wc -l) donors, $(bcftools view -H "$out" | wc -l) sites -> $out"
