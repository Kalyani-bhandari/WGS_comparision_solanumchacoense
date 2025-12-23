#!/bin/bash
set -e

# Load environment
source ~/.bashrc
conda activate wgs_env

# Output directory
OUTDIR="results/fastqc_raw"
mkdir -p "$OUTDIR"

echo "Running FastQC on raw FASTQ files..."

# Read sample table
while read -r sample group R1 R2; do
  # Skip header line
  [[ "$sample" == "sample" ]] && continue

  echo "Processing sample: $sample"
  echo "  R1: $R1"
  echo "  R2: $R2"

  fastqc -t 4 "$R1" "$R2" -o "$OUTDIR"

done < config/samples.tsv

echo "FastQC completed for all samples."

