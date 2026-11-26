#!/bin/bash
# Keep the records of a site VCF that souporcell's covered-variant filter keeps:
# those tested at a position where MIN_COV <= read depth < MAX_COV in BAM.
#
# The output (header, records, record order and repeated records) is that of
# the whole-BAM pipe souporcell runs, which this script replaces. As the rule
# ran it before:
#   samtools depth -d 0 BAM \
#       | awk -v lo=MIN_COV -v hi=MAX_COV '{ if ($3 >= lo && $3 < hi) print $1"\t"$2"\t"($2+1)"\t"$3 }' \
#       | bedtools merge -i - >merged.bed
#   bcftools view -h VCF | grep -v '^##bcftools_'
#   bedtools intersect -wa -a VCF -b merged.bed
# That pipe writes one line for every covered base of the BAM. Here samtools
# reports depth only at the positions the records are tested at
# (`samtools depth -b`), and awk picks the records.
#
# What the old pipe did, and this script reproduces:
#  * The tested positions. The awk step wrote depth position p as the BED
#    interval (p, p + 1), which is base p + 1, and bedtools places a record at
#    POS with a REF of length L on bases POS .. POS+L-1. A record is therefore
#    tested at depth positions POS-1 .. POS+L-2, one base before itself, as in
#    souporcell. A SNV at chr1:1000 is kept when chr1:999 is in range.
#  * Repeated records. bedtools merge joins consecutive in-range positions into
#    one interval, and intersect -wa writes a record once for each interval it
#    overlaps. A deletion at 1000 whose REF has 4 bases is tested at
#    999 .. 1002. With 999 and 1000 in range, 1001 out of range and 1002 in
#    range, it is written twice.
#  * The header, as htslib writes it (htsfile -ch prints what bcftools view -h
#    prints) without the ##bcftools_ lines, which carry a date. htslib adds a
#    ##FILTER=<ID=PASS> line when the file has none.
#  * Uncovered positions. samtools depth (1.9 in the old pipe, 1.15 or 1.22
#    here) does not report a position that no read spans, so such a position
#    is never in range, whatever MIN_COV. All three report a position spanned
#    only by a deletion or a splice gap, at depth 0.
#
# Records this script refuses (the job stops with a message):
#  * a symbolic first ALT (<DEL>, <INS>, ...), which bedtools sized from INFO
#    END= or SVLEN=, whichever came first;
#  * POS < 1, which bedtools refused too.
#
# Known differences from the old pipe, all on inputs the workflow's panels do
# not have:
#  * CRLF line ends: bedtools dropped the \r, this script keeps it.
#  * A blank line in a plain-text VCF: the old pipe skipped it, this script
#    stops (with the POS < 1 message). htslib refuses such a file anyway.
#  * Contigs named browser... or track...: bedtools read those records as
#    header lines and dropped them; this script keeps them.
#  * ##META header lines with a Values=[...] list: bcftools 1.9 dropped them,
#    htsfile keeps them.
#  * A POS with a leading zero: BusyBox awk 1.32 (the samtools 1.15 image) reads
#    it as octal.
#
# Speed: samtools 1.22 took a median of about 8 CPU-s per GB of BAM on
# production captures. samtools 1.15 gives the same output 11 to 38 times more
# slowly on panels of about 200,000 sites and on the test fixture, and 13 per
# cent more slowly on a 6-million-site panel; 1.16 to 1.21 were not timed. The
# log notes a samtools older than 1.22.
#
# Usage: covered_variants.sh BAM VCF MIN_COV MAX_COV OUT_VCF
set -euo pipefail
bam=$1
vcf=$2
min_cov=$3
max_cov=$4
out=$5
bed=$out.tested.bed
pos=$out.in_range.tsv
trap 'rm -f "$bed" "$pos"' EXIT

samtools --version | awk 'NR == 1 {
    print
    split($2, v, ".")
    if (v[1] + 0 == 1 && v[2] + 0 < 22)
        print "note: samtools " $2 " is older than 1.22: depth -b gives the same output, more slowly"
}' >&2

# bgzip -dc refuses plain text, so the reader is chosen from the file's first
# two bytes (1f 8b for gzip and bgzip). Reading from stdin skips the file-name
# suffix check of older bgzip (htslib 1.9 refuses a .bgz file).
if [[ "$(od -An -N2 -tx1 "$vcf" | tr -d ' \n')" == 1f8b ]]; then
    read_vcf() { bgzip -dc <"$vcf"; }
else
    read_vcf() { cat "$vcf"; }
fi

# 1. The depth positions each record is tested at, POS-1 .. POS+L-2, as BED
#    (0-based, half-open).
read_vcf | awk -F'\t' '
    /^#/ { next }
    $2 < 1 || $5 ~ /^</ {
        rec = $1 ":" $2 " " $4 ">" $5
        print "covered_variants.sh: cannot test " rec ": symbolic ALT or POS < 1." > "/dev/stderr"
        exit 1
    }
    {
        s = $2 - 2
        e = s + length($4)
        if (s < 0) s = 0
        if (e > s) print $1 "\t" s "\t" e
    }' >"$bed"

# 2. Depth at those positions only, kept when MIN_COV <= depth < MAX_COV.
#    -d 0 lifts the depth cap of samtools before 1.13; later versions have no
#    cap and ignore the option.
if [[ -s "$bed" ]]; then
    samtools depth -d 0 -b "$bed" "$bam" |
        awk -F'\t' -v lo="$min_cov" -v hi="$max_cov" \
            '$3 + 0 >= lo + 0 && $3 + 0 < hi + 0 { print $1 "\t" $2 }' >"$pos"
else
    : >"$pos"
fi

# 3. The header, then each record once for each run of consecutive in-range
#    positions among its tested positions (each run is one merged interval
#    that bedtools would have matched).
{
    htsfile -ch "$vcf" | grep -v '^##bcftools_'
    read_vcf | awk -F'\t' -v pos="$pos" '
        BEGIN { while ((getline line < pos) > 0) ok[line] = 1 }
        /^#/ { next }
        {
            s = $2 - 1
            e = $2 - 2 + length($4)
            if (s < 1) s = 1
            n = 0
            for (q = s; q <= e; q++)
                if ((($1 "\t" q) in ok) && (q == s || !(($1 "\t" (q - 1)) in ok)))
                    n++
            for (i = 0; i < n; i++) print
        }'
} >"$out"
