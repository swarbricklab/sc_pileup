#!/bin/bash
# Touch empty placeholder inputs so `snakemake -n` can build the full DAG in CI
# without any real data (all three engines, one short capture).
set -eu
base=.tests/ci
touch $base/pileup.bam $base/pileup.bam.bai $base/barcodes.tsv.gz \
      $base/sites.vcf.gz $base/genome.fa $base/genome.fa.fai
echo "placeholders touched"
