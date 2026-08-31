#!/usr/bin/env bash
# build_informative_panel.sh -- construct per-pool demultiplexing SNP panels.
#
# Reusable, project-agnostic tool. It lives with sc_pileup but is NOT part of the
# pileup workflow: it is called by (a) the workflow's optional panel rule for
# standalone use, and (b) orchestration DVC stages (e.g. a bcarc_chromium
# `build_informative_panels` stage). Pure: every project-specific input is an
# explicit argument; it reads no project config.
#
# Given a donor-genotype VCF (cohort or project), it emits, per capture, the SNP
# panel the demux pileup tallies at, applying up to three COMPOSABLE, optional
# steps -- all folded into one `bcftools view`:
#   [--common-sites REF]  restrict to a common-SNP backbone (harmonises WGS vs
#                         microarray-imputed donors onto identical sites). REF is
#                         a bgzipped+indexed VCF or a BED/tab sites file.
#   --donors  d1,d2,...   subset to the pool's donors (REQUIRED).
#   [--informative]       keep only sites polymorphic among those donors
#                         (`bcftools -c 1:minor`: drops all-hom-ref/all-hom-alt).
#                         The allele count is recomputed on the donor subset, so
#                         this is pool-specific.
#
# Purpose presets (set by the caller, not here):
#   demux    : --common-sites <backbone> --informative
#   Track B  : --common-sites <backbone>                 (germline hets, no info filter)
#   Track A  : (neither)                                 (somatic; keep every site)
#
# Modes:
#   single : --capture ID --donors LIST
#   batch  : --captures TSV   (lines: "<capture_id>\t<d1;d2;...>"; '#' comments ok)
# Output : <out-dir>/<capture_id>/informative.vcf.gz (+ .csi)
#
# Assumes the input VCF already carries proper donor GTs (source merge / harmonise
# and missing->0/0 are done upstream). Requires bcftools on PATH (run inside the
# demuxafy / bcftools container).

set -euo pipefail

vcf=""; out_dir=""; common=""; informative=0
captures_tsv=""; capture=""; donors=""

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --vcf)          vcf=$2; shift 2;;
    --out-dir)      out_dir=$2; shift 2;;
    --common-sites) common=$2; shift 2;;
    --informative)  informative=1; shift;;
    --captures)     captures_tsv=$2; shift 2;;
    --capture)      capture=$2; shift 2;;
    --donors)       donors=$2; shift 2;;
    -h|--help)      usage 0;;
    *) echo "unknown arg: $1" >&2; usage 1;;
  esac
done

[ -n "$vcf" ] && [ -n "$out_dir" ] || { echo "need --vcf and --out-dir" >&2; usage 1; }

build_one() {
  # $1=capture_id  $2=donors (',' or ';' separated)
  local cap=$1 dl=${2//;/,}
  local od="$out_dir/$cap"
  mkdir -p "$od"
  local args=(view -s "$dl" --force-samples)
  [ -n "$common" ]        && args+=(-T "$common")   # streaming site intersect
  [ "$informative" = 1 ]  && args+=(-c 1:minor)     # recomputed on the subset
  args+=("$vcf" -Oz -o "$od/informative.vcf.gz")
  bcftools "${args[@]}"
  bcftools index -f "$od/informative.vcf.gz"
  echo "  $cap: $(bcftools view -H "$od/informative.vcf.gz" | wc -l) sites, donors=[$dl]"
}

if [ -n "$captures_tsv" ]; then
  echo "batch panel build from $captures_tsv -> $out_dir"
  while IFS=$'\t' read -r cap dl _; do
    [ -z "${cap:-}" ] && continue
    case "$cap" in \#*) continue;; esac
    build_one "$cap" "$dl"
  done < "$captures_tsv"
elif [ -n "$capture" ] && [ -n "$donors" ]; then
  build_one "$capture" "$donors"
else
  echo "need either --captures TSV or (--capture ID --donors LIST)" >&2; usage 1
fi
