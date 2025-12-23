#!/bin/bash
set -euo pipefail

#Environment
source ~/.bashrc
conda activate wgs_env

#Load_shared_paths
source config/paths.sh

#Directories
BAMDIR="$BAM_DIR"
VCFDIR="$VCF_DIR"

mkdir -p "$VCFDIR"

echo "Starting variant calling using bcftools..."

#Loop_over_samples
while read -r sample group R1 R2; do
  # Skip header
  [[ "$sample" == "sample" ]] && continue

  BAM="$BAMDIR/${sample}.sorted.bam"
  OUTVCF="$VCFDIR/${sample}.vcf.gz"

  echo "--------------------------------------------"
  echo "Sample: $sample"
  echo "Input BAM: $BAM"
  echo "Output VCF: $OUTVCF"

  #Check
  if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found for $sample"
    exit 1
  fi

  #Variant_calling
  bcftools mpileup -Ou -f "$REF" "$BAM" \
    | bcftools call -mv -Oz -o "$OUTVCF"

  #Index_VCF
  bcftools index "$OUTVCF"

done < config/samples.tsv

echo "Variant calling completed for all samples."

