#!/bin/bash
set -euo pipefail

source ~/.bashrc
conda activate wgs_env

ANNOTDIR="results/snpeff_annotated"
OUTDIR="results/snpeff_filtered"

mkdir -p "$OUTDIR"

echo "Extracting HIGH and MODERATE impact variants..."

for GROUP in Resistance_specific Susceptible_specific; do
  VCF="$ANNOTDIR/${GROUP}.annotated.vcf.gz"

  # HIGH impact
  bcftools view -i 'INFO/ANN ~ "HIGH"' "$VCF" -Oz \
    -o "$OUTDIR/${GROUP}.HIGH.vcf.gz"

  bcftools index "$OUTDIR/${GROUP}.HIGH.vcf.gz"

  # MODERATE impact
  bcftools view -i 'INFO/ANN ~ "MODERATE"' "$VCF" -Oz \
    -o "$OUTDIR/${GROUP}.MODERATE.vcf.gz"

  bcftools index "$OUTDIR/${GROUP}.MODERATE.vcf.gz"

done

echo "Impact-based filtering completed."

