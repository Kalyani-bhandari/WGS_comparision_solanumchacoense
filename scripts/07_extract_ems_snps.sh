#!/bin/bash
set -euo pipefail

# ================================
# Environment
# ================================
source ~/.bashrc
conda activate wgs_env

# ================================
# Directories
# ================================
FILTERDIR="results/vcf_filtered"
OUTDIR="results/EMS_vcfs"

mkdir -p "$OUTDIR"

echo "Extracting EMS-specific SNPs (G>A, C>T)..."

# ================================
# Loop over filtered VCFs
# ================================
for VCF in "$FILTERDIR"/*.filtered.vcf.gz; do
  BASENAME=$(basename "$VCF" .filtered.vcf.gz)
  OUTVCF="$OUTDIR/${BASENAME}.EMS.vcf.gz"

  echo "--------------------------------------------"
  echo "Processing: $BASENAME"

  bcftools view \
    -i 'TYPE="snp" && ((REF="G" && ALT="A") || (REF="C" && ALT="T"))' \
    "$VCF" -Oz -o "$OUTVCF"

  bcftools index "$OUTVCF"

done

echo "EMS SNP extraction completed."

