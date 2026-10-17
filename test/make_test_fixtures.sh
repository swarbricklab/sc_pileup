#! /bin/bash
# Regenerate the sc_pileup test fixtures (test/data/) from PUBLIC sources.
#
# The fixtures are DVC-tracked; this script documents their provenance and lets
# anyone with access to the NCI a56 inputs rebuild them deterministically. It is
# NOT run by CI — CI consumes the DVC-tracked outputs.
#
# Provenance
# ----------
#  * BAM + barcodes: a tiny_pbmc capture from swarbricklab/chromium-preprocessing's
#    PUBLIC test data (downsampled 10x PBMC, run through Cell Ranger count). We take
#    one capture (PBMC_001) and subset its cellranger BAM to two windows:
#        chr19:1-10,000,000  and  chr22:20,000,000-30,000,000
#    (gene-dense, well-covered; two chromosomes so popscle's per-chromosome
#    split/merge is exercised).
#  * sites VCF: the souporcell IGSR/1000G GRCh38 common-SNP panel
#    (common_variants_grch38.vcf), subset to the same two windows. The panel's data
#    rows use Ensembl-style contig names (1,19,22) while the cellranger BAM uses
#    UCSC-style (chr1,chr19); we rename the kept rows to chr19/chr22 to match.
#    The pileup is genotype-FREE, so the VCF only needs to be a set of positions.
#
# Output: test/data/{pileup.bam,pileup.bam.bai,barcodes.tsv.gz,sites.vcf.gz,sites.vcf.gz.tbi}
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate snakemake_7.32.4
module load samtools bcftools htslib 2>/dev/null || true

# -- Public NCI-resident inputs (override via env if your paths differ) --
SRC_BAM=${SRC_BAM:-/scratch/a56/jr9959/chromium-preprocessing-pubprep/test/results/cellranger/count_test/bam/count/captures/PBMC_001/assigned_alignments.bam}
SRC_BC=${SRC_BC:-/scratch/a56/jr9959/chromium-preprocessing-pubprep/test/results/cellranger/count_test/count/filtered/PBMC_001/barcodes.tsv.gz}
SRC_VCF=${SRC_VCF:-/g/data/a56/backups/sunwu/souporcell/genotyping/annotation/common_variants_grch38.vcf}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/data"
mkdir -p "$OUT"

echo ">> subset BAM (chr19:1-10M + chr22:20-30M)"
samtools view -b -o "$OUT/pileup.bam" "$SRC_BAM" chr19:1-10000000 chr22:20000000-30000000
samtools index "$OUT/pileup.bam"

echo ">> barcodes (filtered cell list, as-is)"
install -m 0644 "$SRC_BC" "$OUT/barcodes.tsv.gz"

echo ">> sites VCF: two pileup windows + >=1 padding site on EVERY other autosome"
# popscle dsc-pileup runs once per autosome chr1..chr22 and FATALs on any
# chromosome whose (split) VCF has zero variants. Production panels are
# genome-wide so this never bites; our two-window fixture would leave 20 empty
# chromosomes. So we keep the real content on chr19/chr22 and add a few padding
# sites (PAD per chromosome) on every other autosome. The BAM still only has
# reads on chr19/chr22 — dsc-pileup tolerates sites-with-zero-reads fine.
PAD=${PAD:-50}
HDR=$(mktemp); REC=$(mktemp)
# records (chrom renamed to chr-style), from the karyotypic-named source
awk -F'\t' -v pad="$PAD" 'BEGIN{OFS="\t"} /^#/{next}
     $1 !~ /^[0-9]+$/ {next}                               # autosomes only
     $1=="19" && $2<=10000000              {$1="chr19"; print; next}
     $1=="22" && $2>=20000000 && $2<=30000000 {$1="chr22"; print; next}
     n[$1] < pad                           {c=$1; $1="chr"c; print; n[c]++}' "$SRC_VCF" > "$REC"
# header: source meta (minus contig) + ##contig in BAM @SQ order + #CHROM
grep '^##' "$SRC_VCF" | grep -v '^##contig' > "$HDR"
samtools view -H "$OUT/pileup.bam" | awk -F'\t' '/^@SQ/{name="";len="";
    for(i=1;i<=NF;i++){if($i~/^SN:/)name=substr($i,4); if($i~/^LN:/)len=substr($i,4)}
    print "##contig=<ID="name",length="len">"}' >> "$HDR"
grep -m1 '^#CHROM' "$SRC_VCF" >> "$HDR"
# bcftools sort orders records by the header's ##contig order (= BAM @SQ order),
# giving a VCF whose contig ordering matches the BAM (popscle's requirement).
cat "$HDR" "$REC" | bcftools sort -Oz -o "$OUT/sites.vcf.gz" 2>/dev/null
rm -f "$HDR" "$REC"
tabix -p vcf "$OUT/sites.vcf.gz"

echo ">> done:"; du -h "$OUT"/* | sort -h
echo ">> remember to: dvc add test/data && dvc push"
