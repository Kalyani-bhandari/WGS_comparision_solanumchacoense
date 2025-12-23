#!/bin/bash
set -euo pipefail

# ================================
# Environment
# ================================
source ~/.bashrc
conda activate wgs_env

# ================================
# Load paths
# ================================
source config/paths.sh

# ================================
# Directories
# ================================
VCFDIR="$VCF_DIR"
OUTDIR="results/vcf_filtered"

mkdir -p "$OUTDIR"

echo "Starting VCF quality filtering..."

# ================================
# Loop over samples
# ================================
while read -r sample group R1 R2; do
  [[ "$sample" == "sample" ]] && continue

  INVCF="$VCFDIR/${sample}.vcf.gz"
  OUTVCF="$OUTDIR/${sample}.filtered.vcf.gz"

  echo "--------------------------------------------"
  echo "Filtering sample: $sample"
  echo "Input:  $INVCF"
  echo "Output: $OUTVCF"

  if [[ ! -f "$INVCF" ]]; then
    echo "ERROR: VCF not found for $sample"
    exit 1
  fi

  bcftools filter \
    -e 'QUAL<30 || DP<8 || MQ<30 || MQ0F>0.05 || VDB<0.0001 || SGB<-0.6' \
    "$INVCF" -Oz -o "$OUTVCF"

  bcftools index "$OUTVCF"

done < config/samples.tsv

echo "VCF filtering completed."

