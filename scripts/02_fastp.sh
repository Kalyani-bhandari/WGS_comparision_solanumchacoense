#!/bin/bash
set -e

#Load_environment
source ~/.bashrc
conda activate wgs_env

#Output_directory
OUTDIR="results/fastp"
mkdir -p "$OUTDIR"

echo "Running fastp on all samples..."

#Read_sample_table
while read -r sample group R1 R2; do
  # Skip header
  [[ "$sample" == "sample" ]] && continue

  echo "Processing sample: $sample"
  echo "  R1: $R1"
  echo "  R2: $R2"

  fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$OUTDIR/${sample}_R1.trimmed.fastq.gz" \
    -O "$OUTDIR/${sample}_R2.trimmed.fastq.gz" \
    --html "$OUTDIR/${sample}.html" \
    --json "$OUTDIR/${sample}.json" \
    --detect_adapter_for_pe \
    --thread 8

done < config/samples.tsv

echo "fastp trimming completed for all samples."

