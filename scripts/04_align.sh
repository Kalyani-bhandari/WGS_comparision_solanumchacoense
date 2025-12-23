#!/bin/bash
set -e

#Load_environment
source ~/.bashrc
conda activate wgs_env

#Load_paths
source config/paths.sh

#Input_output_directories
TRIMDIR="results/fastp"
BAMDIR="results/bam"

mkdir -p "$BAMDIR"

echo "Starting BWA alignment for all samples..."

#Read_sample_table
while read -r sample group R1 R2; do

  #Skip_header
  [[ "$sample" == "sample" ]] && continue

  R1_TRIM="$TRIMDIR/${sample}_R1.trimmed.fastq.gz"
  R2_TRIM="$TRIMDIR/${sample}_R2.trimmed.fastq.gz"
  OUTBAM="$BAMDIR/${sample}.sorted.bam"

  echo "----------------------------------------"
  echo "Aligning sample: $sample"
  echo "R1: $R1_TRIM"
  echo "R2: $R2_TRIM"
  echo "Output: $OUTBAM"

  bwa mem -t 8 "$REF" "$R1_TRIM" "$R2_TRIM" \
    | samtools view -b -@ 4 - \
    | samtools sort -@ 4 -o "$OUTBAM" -

  samtools index "$OUTBAM"

done < config/samples.tsv

echo "Alignment completed for all samples."

