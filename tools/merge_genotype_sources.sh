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
#   [--fill-ref]          set each source's OWN missing GTs to 0/0 BEFORE merging
#                         (bcftools +setGT; a no-call within a source = hom-ref,
#                         valid for WGS). Cross-source missing (a donor's source
#                         lacks a site) stays ./. -- it is NOT invented as hom-ref
#                         (which would falsely assert e.g. a microarray donor is
#                         hom-ref at sites the array never measured).
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
  # --fill-ref is applied PER SOURCE (a no-call within a source = hom-ref, valid
  # for WGS). Cross-source missing (a donor's source lacks a site entirely) is
  # left ./. after the merge -- NOT invented as hom-ref, which would falsely
  # assert e.g. a microarray donor is hom-ref at sites the array never measured.
  if [ "$fill" = 1 ]; then
    bcftools "${args[@]}" "$vcf" -Ou | bcftools +setGT -Oz -o "$p" -- -t . -n 0
  else
    bcftools "${args[@]}" "$vcf" -Oz -o "$p"
  fi
  bcftools index "$p"
  parts+=("$p"); i=$((i+1))
done

if [ ${#parts[@]} -eq 1 ]; then
  bcftools view "${parts[0]}" -Oz -o "$out"
else
  bcftools merge "${parts[@]}" -Oz -o "$out"
fi
bcftools index -f "$out"
echo "merged $(bcftools query -l "$out" | wc -l) donors, $(bcftools view -H "$out" | wc -l) sites -> $out"
