#!/bin/bash
# vartrix at the somatic panel on the ORIGINAL (skip_remap) vs REMAPPED bam.
# Same tool, same sites/barcodes -> matrices align element-wise. Run in
# cumulusprod/souporcell:2.5 (Demuxafy vartrix panics).
set -euo pipefail
SIF=/scratch/a56/jr9959/sc_pileup/resources/souporcell_2.5.sif
FA=/scratch/a56/jr9959/bcarc_chromium/resources/transcriptomes/refdata-gex-GRCh38-2020-A/fasta/genome.fa
BC_GZ=/scratch/a56/jr9959/bcarc_chromium/data/cellranger/count/filtered/Neomet_Pool028a/barcodes.tsv.gz
ORIG=/scratch/a56/jr9959/bcarc_chromium/data/cellranger/bam/count/captures/Neomet_Pool028a/assigned_alignments.bam
REMAP=/scratch/a56/jr9959/sc_pileup/results/souporcell/short/Neomet_Pool028a/remap/minimap_tagged_sorted.bam
SOMVCF=/scratch/a56/jr9959/sc_pileup/resources/somatic/merged_somatic.vcf.gz
OUT=/scratch/a56/jr9959/sc_pileup/results/vartrix_compare/Neomet_Pool028a
THREADS="${1:-16}"

mkdir -p "$OUT"
# Restrict the panel to contigs present in the cellranger reference FASTA —
# DRAGEN somatic VCFs call on alt/decoy contigs (chrUn_*, *_alt, *_random) that
# refdata-gex-GRCh38-2020-A lacks; vartrix errors on those, and they never have
# scRNA coverage anyway. (Panel build in somatifind should do the same.)
zcat "$SOMVCF" | awk -v fai="$FA.fai" '
  BEGIN{ while((getline l < fai)>0){ split(l,a,"\t"); keep[a[1]]=1 } }
  /^#/{print; next} ($1 in keep){print}' > "$OUT/somatic.vcf"
zcat "$BC_GZ"  > "$OUT/barcodes.tsv"
echo "panel sites after contig filter: $(grep -vc '^#' "$OUT/somatic.vcf")"

module load singularity 2>/dev/null || true
run_vartrix () {  # $1=bam  $2=tag
  echo "=== vartrix on $2 ($1) ==="
  singularity exec --cleanenv --bind /scratch/a56,/g/data/a56 "$SIF" \
    vartrix --mapq 30 -b "$1" -c "$OUT/barcodes.tsv" --scoring-method coverage \
      --threads "$THREADS" --ref-matrix "$OUT/ref_$2.mtx" --out-matrix "$OUT/alt_$2.mtx" \
      -v "$OUT/somatic.vcf" --fasta "$FA" --umi
}
run_vartrix "$ORIG"  orig
run_vartrix "$REMAP" remap
echo "DONE"; ls -lh "$OUT"/*.mtx | awk '{print $5, $NF}'
