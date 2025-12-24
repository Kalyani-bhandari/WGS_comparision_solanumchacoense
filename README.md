# Whole Genome Sequencing Analysis of EMS Induced Mutations in *Solanum chacoense*

## Overview
This repository contains a whole genome sequencing analysis pipeline developed to identify **EMS induced SNPs** associated with **soft rot resistance and susceptibility** in *Solanum chacoense* M6 mutant lines.

The pipeline compares resistant and susceptible EMS mutant lines to identify:
- EMS mutations (G-A, C-T)
- Resistance specific and susceptibility specific variants
- Shared and unique SNPs across independent mutant lines
- Functionally annotated candidate mutations linked to resistance associated genes

This work supports downstream genetic and functional analyses of Illumina sequencing data.

---

## Pipeline Summary
The analysis workflow includes the following steps:

1. Quality control of raw paired-end reads (FastQC)
2. Adapter trimming and quality filtering (fastp)
3. Reference-guided alignment to the *S. chacoense* M6 genome (BWA-MEM)
4. BAM processing and indexing (SAMtools)
5. Variant calling (bcftools)
6. Variant quality filtering
7. Extraction of EMS-specific SNPs (G-A, C-T)
8. Comparative analysis between resistant and susceptible lines
9. Identification of shared and line-specific SNPs
10. Functional annotation of candidate variants (SnpEff / bcftools csq)

---

## Input Data
- Paired end FASTQ files from EMS mutant lines
- *Solanum chacoense* M6 reference genome (`reference.fa`)
- Gene annotation file (`.gff3`)
- Sample metadata file describing resistant and susceptible lines

> **Note:** Raw sequencing data are not included in this repository.

---

## Software and Tools
The pipeline uses the following bioinformatics tools:

- FastQC
- fastp
- BWA
- SAMtools
- BCFtools
- SnpEff

All tools are intended to be run in a Conda environment.

---

## Repository Structure (under development)
├── config/ # Configuration files (paths, samples)
├── scripts/ # Modular analysis scripts
├── workflow/ # Master pipeline launcher
├── results/ # Output directories (ignored by git)
├── docs/ # Additional documentation

---

## Reproducibility
This repository is being actively refactored into a modular, loop-based, and reproducible pipeline.  
Hard-coded paths and sample names are being replaced with configuration files to improve scalability and portability across computing environments (HPC and local systems).

---

## Author
**Kalyani Bhandari**  
M.S. Student, Plant Pathology  
Montana State University

---

## License
This repository is intended for academic and research use.

