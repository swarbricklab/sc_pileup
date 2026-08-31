#!/bin/bash
# Ground-truth: single unsplit popscle dsc-pileup on the WHOLE BAM/VCF (single
# threaded), equivalent params to the existing chrom-split+merge result for
# Neomet_Pool028a. Same popscle binary (Demuxafy.sif), same BAM, same panel
# (restricted to chr1-22, which is all the split-merge covered), same barcodes,
# default dsc-pileup params (min-MQ 20, cap-BQ 40, excl-flag 3844, CB/UB).
set -euo pipefail
SIF=/scratch/a56/jr9959/bcarc_chromium/resources/Demuxafy.sif
BAM=/scratch/a56/jr9959/bcarc_chromium/data/cellranger/bam/count/captures/Neomet_Pool028a/assigned_alignments.bam
VCF=/scratch/a56/jr9959/bcarc_chromium/data/snp/vcf/donor_keyed.vcf.gz
BC_GZ=/scratch/a56/jr9959/bcarc_chromium/data/cellranger/count/filtered/Neomet_Pool028a/barcodes.tsv.gz
OUT=/scratch/a56/jr9959/sc_pileup/results/popscle_unsplit/Neomet_Pool028a
mkdir -p "$OUT"
CHRS=$(printf 'chr%s,' {1..22}); CHRS=${CHRS%,}
module load singularity 2>/dev/null || true
singularity exec --cleanenv --bind /scratch/a56,/g/data/a56 "$SIF" bash -c "
  set -euo pipefail
  bcftools view -t $CHRS '$VCF' -Oz -o '$OUT/panel_chr1_22.vcf.gz'
  bcftools index '$OUT/panel_chr1_22.vcf.gz'
  zcat '$BC_GZ' > '$OUT/barcodes.tsv'
  echo \"panel sites (chr1-22): \$(zcat '$OUT/panel_chr1_22.vcf.gz' | grep -vc '^#')\"
  echo \"barcodes: \$(wc -l < '$OUT/barcodes.tsv')\"
  echo 'starting unsplit dsc-pileup (single-threaded, may take many hours)...'
  time popscle dsc-pileup --sam '$BAM' --vcf '$OUT/panel_chr1_22.vcf.gz' \
       --group-list '$OUT/barcodes.tsv' --out '$OUT/unsplit'
"
echo DONE; ls -lh "$OUT"/unsplit.* 2>/dev/null
