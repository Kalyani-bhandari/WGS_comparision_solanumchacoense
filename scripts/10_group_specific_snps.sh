#!/bin/bash
set -euo pipefail

source ~/.bashrc
conda activate wgs_env

OUTDIR="results/group_specific"
SHAREDDIR="results/shared_groups"

mkdir -p "$OUTDIR"

echo "Identifying group-specific EMS SNPs..."

# Resistance-specific
bcftools isec \
  "$SHAREDDIR/Resistant_shared/0000.vcf.gz" \
  "$SHAREDDIR/Susceptible_shared/0000.vcf.gz" \
  -p "$OUTDIR/Resistance_specific"

# Susceptible-specific
bcftools isec \
  "$SHAREDDIR/Susceptible_shared/0000.vcf.gz" \
  "$SHAREDDIR/Resistant_shared/0000.vcf.gz" \
  -p "$OUTDIR/Susceptible_specific"

# Compress + index
for dir in "$OUTDIR"/*; do
  for f in "$dir"/*.vcf; do
    bgzip "$f"
    bcftools index "$f.gz"
  done
done

echo "Group-specific SNP extraction completed."

