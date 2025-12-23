#!/bin/bash
set -euo pipefail

source ~/.bashrc
conda activate wgs_env
source config/paths.sh

INBASE="results/group_specific"
OUTBASE="results/snpeff_annotated"

mkdir -p "$OUTBASE"

echo "Annotating EMS SNPs with SnpEff..."

for GROUP in Resistance_specific Susceptible_specific; do
  INVCF="$INBASE/$GROUP/0000.vcf.gz"
  OUTVCF="$OUTBASE/${GROUP}.annotated.vcf.gz"

  echo "--------------------------------------------"
  echo "Annotating: $GROUP"

  snpEff -noStats solanum_chacoense_m6 "$INVCF" \
    | bgzip > "$OUTVCF"

  bcftools index "$OUTVCF"
done

echo "SnpEff annotation completed."

