#!/bin/bash

# ================================
# Reference genome and annotation
# ================================

# Update these paths on the HPC system (Tempest)
REF="/home/g89x126/sc_wgs/reference.fa"
GFF="/home/g89x126/sc_wgs/m6.hc.gene_models.all_chrs_corrected.gff3"

# ================================
# Output directories (relative to repo root)
# ================================

FASTP_DIR="results/fastp"
FASTQC_RAW_DIR="results/fastqc_raw"
FASTQC_TRIM_DIR="results/fastqc_trimmed"
BAM_DIR="results/bam"
VCF_DIR="results/vcf"

