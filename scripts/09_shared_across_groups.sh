#!/bin/bash
set -euo pipefail

source ~/.bashrc
conda activate wgs_env

EMSDIR="results/EMS_vcfs"
OUTDIR="results/shared_groups"

mkdir -p "$OUTDIR"

echo "Identifying shared EMS SNPs across groups..."

# ================================
# Resistant
# ================================
REST_VCFS=$(ls $EMSDIR/*-Rest.EMS.vcf.gz)

bcftools isec -n=2 $REST_VCFS -p "$OUTDIR/Resistant_shared"

for f in "$OUTDIR/Resistant_shared"/*.vcf; do
  bgzip "$f"
  bcftools index "$f.gz"
done

# ================================
# Susceptible
# ================================
SUSP_VCFS=$(ls $EMSDIR/*-Susp.EMS.vcf.gz)

bcftools isec -n=2 $SUSP_VCFS -p "$OUTDIR/Susceptible_shared"

for f in "$OUTDIR/Susceptible_shared"/*.vcf; do
  bgzip "$f"
  bcftools index "$f.gz"
done

echo "Shared group comparisons completed."

