#!/bin/bash
#SBATCH --job-name=fastp_qc
#SBATCH --output=fastp_qc_%A.out
#SBATCH --error=fastp_qc_%A.err
#SBATCH --time=12:00:00
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate wgs_env

module load fastp

RAW_DIR="/home/g89x126/sc_wgs"
OUT_DIR="${RAW_DIR}/fastp_out"
PNG_DIR="${RAW_DIR}/fastp_png"

mkdir -p $OUT_DIR $PNG_DIR

# Loop through R1 files
for R1 in ${RAW_DIR}/*_R1_001.fastq; do
    SAMPLE=$(basename $R1 _R1_001.fastq)
    R2="${RAW_DIR}/${SAMPLE}_R2_001.fastq"

    echo "Running fastp for: $SAMPLE"

    fastp \
        -i $R1 \
        -I $R2 \
        -o ${OUT_DIR}/${SAMPLE}_R1.trimmed.fastq \
        -O ${OUT_DIR}/${SAMPLE}_R2.trimmed.fastq \
        -h ${OUT_DIR}/${SAMPLE}.html \
        -j ${OUT_DIR}/${SAMPLE}.json \
        --detect_adapter_for_pe \
        --thread 4

    echo "Extracting PNG images for $SAMPLE"
    
    # Extract PNGs from HTML
    mkdir -p "${PNG_DIR}/${SAMPLE}"

    python3 - <<EOF
import os, base64, re

html = open("${OUT_DIR}/${SAMPLE}.html").read()
outdir = "${PNG_DIR}/${SAMPLE}"

# Find all base64 PNG images inside HTML
matches = re.findall(r'data:image/png;base64,([A-Za-z0-9+/=]+)', html)

for i, img in enumerate(matches):
    img_bytes = base64.b64decode(img)
    with open(f"{outdir}/{SAMPLE}_img_{i+1}.png", "wb") as f:
        f.write(img_bytes)
EOF

done
