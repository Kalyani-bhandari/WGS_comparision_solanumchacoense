#!/bin/bash
set -euo pipefail

echo " WGS EMS Mutation Discovery Pipeline"

#QC_and_preprocessing

echo "Step 1: FastQC on raw reads"
bash scripts/01_fastqc.sh

echo "Step 2: Read trimming with fastp"
bash scripts/02_fastp.sh

echo "Step 3: FastQC on trimmed reads"
bash scripts/03_fastqc_trimmed.sh

#Alignment_and_variant_calling
echo "Step 4: Alignment to reference genome"
bash scripts/04_align.sh

echo "Step 5: Variant calling (bcftools)"
bash scripts/05_variant_calling.sh

#Variant_filtering_and_EMS_extraction
echo "Step 6: Variant quality filtering"
bash scripts/06_filter_vcf.sh

echo "Step 7: EMS SNP extraction (G>A, C>T)"
bash scripts/07_extract_ems_snps.sh

#Comparative_analyses
echo "Step 8: Pairwise comparisons within EMS lines"
bash scripts/08_pairwise_within_lines.sh

echo "Step 9: Shared EMS SNPs across groups"
bash scripts/09_shared_across_groups.sh

echo "Step 10: Group-specific EMS SNPs"
bash scripts/10_group_specific_snps.sh

#Functional_annotation
echo "Step 11: Build SnpEff database (run once if needed)"
bash scripts/11_snpeff_build.sh

echo "Step 12: SnpEff annotation"
bash scripts/12_snpeff_annotate.sh

echo "Step 13: Extract HIGH/MODERATE impact variants"
bash scripts/13_snpeff_extract_impacts.sh

echo " Pipeline completed successfully"

