#!/bin/bash
# Classify the records of a site VCF and report how many are NOT biallelic SNVs.
#
# cellSNP-lite piles up biallelic SNVs, so indels, MNVs, multiallelic and symbolic
# records in a somatic (Track A) site VCF can never be observed; left in, they
# only inflate the apparent "unreachable" fraction downstream. This writes a
# one-row TSV of counts and, when OUT_VCF is given, a bgzipped+indexed VCF
# holding only the biallelic SNVs. MODE (warn|drop) sets the log message and the
# report's non_snv_action column.
#
# Classes (first match wins, REF/ALT compared case-insensitively):
#   multiallelic   ALT contains ','
#   biallelic_snv  REF and ALT are each a single A/C/G/T base
#   mnv            REF and ALT are equal-length (>1) base strings
#   indel          REF and ALT are base strings of different lengths
#   other          anything else (symbolic <DEL>, '*', '.', N, breakends, ...)
#
# Usage: check_site_types.sh VCF MODE(warn|drop) CAPTURE_ID REPORT_TSV [OUT_VCF]
set -euo pipefail
vcf=$1
mode=$2
capture=$3
report=$4
out_vcf=${5:-}

counts=$(
    bcftools view -H "$vcf" | awk -F'\t' '
        BEGIN { n = snv = multi = mnv = indel = other = 0 }
        {
            n++
            ref = toupper($4); alt = toupper($5)
            if (alt ~ /,/) { multi++; next }
            if (ref ~ /^[ACGT]$/ && alt ~ /^[ACGT]$/) { snv++; next }
            if (ref ~ /^[ACGTN]+$/ && alt ~ /^[ACGTN]+$/) {
                if (length(ref) == length(alt)) mnv++; else indel++
                next
            }
            other++
        }
        END { printf "%d\t%d\t%d\t%d\t%d\t%d", n, snv, multi, mnv, indel, other }'
)
read -r n snv multi mnv indel other <<<"$counts"
non_snv=$((n - snv))
if [[ "$mode" == "drop" ]]; then action=dropped; else action=kept; fi

printf 'capture_id\tn_records\tn_biallelic_snv\tn_non_snv\tn_multiallelic\tn_mnv\tn_indel\tn_other\tnon_snv_action\n' >"$report"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$capture" "$n" "$snv" "$non_snv" "$multi" "$mnv" "$indel" "$other" "$action" >>"$report"

echo "site types for $capture ($vcf): $n records, $snv biallelic SNVs," \
    "$non_snv not (multiallelic=$multi mnv=$mnv indel=$indel other=$other)"
if ((non_snv > 0)); then
    if [[ "$mode" == "drop" ]]; then
        echo "dropping $non_snv non-biallelic-SNV record(s) before piling up" \
            "(somatic.non_snv: drop)"
    else
        echo "WARNING: $non_snv of $n record(s) in the somatic site VCF are not" \
            "biallelic SNVs; cellSNP-lite cannot pile them up, so they will" \
            "inflate the unreachable fraction. Set somatic.non_snv: drop to remove" \
            "them before piling up."
    fi
fi

if [[ -n "$out_vcf" ]]; then
    # strip bcftools' own `##bcftools_*` header lines (they carry a Date) so the
    # subset VCF is reproducible; records are filtered, not rewritten.
    (
        bcftools view -h "$vcf" | grep -v '^##bcftools_'
        bcftools view -H "$vcf" | awk -F'\t' \
            'toupper($4) ~ /^[ACGT]$/ && toupper($5) ~ /^[ACGT]$/'
    ) | bcftools view --no-version -Oz -o "$out_vcf" -
    bcftools index -f "$out_vcf"
fi
