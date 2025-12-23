#!/bin/bash
set -e

#Load_environment
source ~/.bashrc
conda activate wgs_env

#Input_and_output_directories
TRIMDIR="results/fastp"
OUTDIR="results/fastqc_trimmed"

mkdir -p "$OUTDIR"

echo "Running FastQC on trimmed FASTQ files..."

#Read_sample_table
while read -r sample group R1 R2; do

  #Skip_header
  [[ "$sample" == "sample" ]] && continue

  R1_TRIM="$TRIMDIR/${sample}_R1.trimmed.fastq.gz"
  R2_TRIM="$TRIMDIR/${sample}_R2.trimmed.fastq.gz"

  echo "Processing sample: $sample"
  echo "  R1 trimmed: $R1_TRIM"
  echo "  R2 trimmed: $R2_TRIM"

  fastqc -t 8 "$R1_TRIM" "$R2_TRIM" -o "$OUTDIR"

done < config/samples.tsv

echo "FastQC on trimmed reads completed."

