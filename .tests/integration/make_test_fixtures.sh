#! /bin/bash
# Regenerate the sc_pileup test fixtures (.tests/integration/data/) from a PUBLIC source.
#
# The fixtures are DVC-tracked; this script documents their provenance and lets
# anyone rebuild them deterministically. It is NOT run by CI — CI consumes the
# DVC-tracked outputs.
#
# Provenance
# ----------
# Demuxafy's public test dataset "TestData4PipelineSmall" — the OneK1K 14-donor
# PBMC pool (Neavin et al., Powell Lab), the canonical multi-donor demultiplexing
# example. Downloaded from the (updated) Google Drive links in the Demuxafy docs
# (https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DataPrep.html):
#
#   wget -O TestData4PipelineSmall.tar.gz \
#     'https://drive.usercontent.google.com/download?id=1cx2Wv46hoYwn659ixSFcYGVKKQw8EXDy&export=download&confirm=t'
#   # md5: d2d5f1c2d405588de8e88b2419eebffa
#
# The tarball carries test_dataset/outs/pooled.sorted.bam (+.bai), the filtered
# barcodes, and test_dataset.vcf (14 donor genotypes, 229,719 variants). The pool
# is 1400 cells; the data are GRCh38p10 with ENSEMBL-style contigs (1,2,...).
#
# The SAME downsampled pool is shared with snp_demux: sc_pileup produces the
# per-cell pileup here, and snp_demux consumes that pileup + the 14-donor genotypes
# (test_dataset.vcf) as the demultiplexing panel.
#
# Transforms
# ----------
#  * BAM: subsample to 10% of reads (seed 42) to keep the fixture CI-small while
#    retaining genome-wide SNP breadth (needed to separate 14 donors), then rename
#    contigs ENSEMBL->UCSC (1->chr1 ... MT->chrM) to match the cellranger-style
#    genome FASTA used by souporcell. Primary-chromosome coordinates are identical
#    across GRCh38 builds, so the rename is coordinate-safe.
#  * sites VCF: the 14-donor VCF, renamed to chr-style, restricted to autosomes
#    (chr1..chr22) and sorted into the BAM's @SQ order (popscle requires the VCF
#    contig order to match the BAM). It KEEPS the 14 donor genotypes (it is the
#    donor panel): cellsnp-lite and popscle use only the positions, but the
#    genotypes flow through souporcell's covered-variants -> variants.vcf, which
#    snp_demux consumes as `--known_genotypes` for genotype-aware clustering. (For
#    a standalone, non-demux pileup a sites-only common-SNP panel would do.)
#
# Output: .tests/integration/data/{pileup.bam,pileup.bam.bai,barcodes.tsv.gz,sites.vcf.gz,sites.vcf.gz.tbi}
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate snakemake_7.32.4
module load samtools bcftools htslib 2>/dev/null || true

# -- Source tarball (download first per the header; override to point elsewhere) --
SRC=${SRC:-../demuxafy_testdata/extracted/TestData4PipelineSmall}
SRC_BAM="$SRC/test_dataset/outs/pooled.sorted.bam"
SRC_BC="$SRC/test_dataset/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/barcodes.tsv"
SRC_VCF="$SRC/test_dataset.vcf"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/data"
mkdir -p "$OUT"

# Ensembl -> UCSC chr-style rename map (primary assembly; coord-identical on GRCh38)
MAP=$(mktemp)
{ for i in $(seq 1 22); do echo "$i chr$i"; done; echo "X chrX"; echo "Y chrY"; echo "MT chrM"; } > "$MAP"

echo ">> BAM: 10% subsample (seed 42) + Ensembl->chr reheader"
samtools view -b -s 42.1 -o "$OUT/_ds.bam" "$SRC_BAM"
samtools view -H "$OUT/_ds.bam" | awk 'BEGIN{OFS="\t"} /^@SQ/{for(i=1;i<=NF;i++){if($i~/^SN:/){n=substr($i,4); if(n~/^([0-9]+|X|Y)$/)$i="SN:chr"n; else if(n=="MT")$i="SN:chrM"}}} {print}' > "$OUT/_hdr.sam"
samtools reheader "$OUT/_hdr.sam" "$OUT/_ds.bam" > "$OUT/pileup.bam"
samtools index "$OUT/pileup.bam"
rm -f "$OUT/_ds.bam" "$OUT/_hdr.sam"

echo ">> barcodes (filtered cell list, as-is)"
gzip -c "$SRC_BC" > "$OUT/barcodes.tsv.gz"

echo ">> sites VCF: 14-donor panel (GT kept), autosomes, chr-style, sorted to BAM @SQ order"
bgzip -c "$SRC_VCF" > "$OUT/_orig.vcf.gz"; tabix -f -p vcf "$OUT/_orig.vcf.gz"
bcftools annotate --rename-chrs "$MAP" "$OUT/_orig.vcf.gz" -Oz -o "$OUT/_chr.vcf.gz"
tabix -f -p vcf "$OUT/_chr.vcf.gz"
REG=$(for i in $(seq 1 22); do printf "chr$i,"; done | sed 's/,$//')
bcftools view -r "$REG" "$OUT/_chr.vcf.gz" -Ov -o "$OUT/_sites.vcf"
grep '^##' "$OUT/_sites.vcf" | grep -v '^##contig' > "$OUT/_hdr"
samtools view -H "$OUT/pileup.bam" | awk '/^@SQ/{name="";len="";
    for(i=1;i<=NF;i++){if($i~/^SN:/)name=substr($i,4); if($i~/^LN:/)len=substr($i,4)}
    if(name!="")print "##contig=<ID="name",length="len">"}' >> "$OUT/_hdr"
grep -m1 '^#CHROM' "$OUT/_sites.vcf" >> "$OUT/_hdr"
grep -v '^#' "$OUT/_sites.vcf" > "$OUT/_recs"
cat "$OUT/_hdr" "$OUT/_recs" | bcftools sort -Oz -o "$OUT/sites.vcf.gz" 2>/dev/null
tabix -f -p vcf "$OUT/sites.vcf.gz"
rm -f "$OUT"/_orig.vcf.gz* "$OUT"/_chr.vcf.gz* "$OUT/_sites.vcf" "$OUT/_hdr" "$OUT/_recs" "$MAP"

echo ">> done:"; du -h "$OUT"/* | sort -h
echo ">> remember to: dvc add .tests/integration/data && dvc push"
