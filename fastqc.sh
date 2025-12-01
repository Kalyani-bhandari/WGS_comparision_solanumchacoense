source ~/.bashrc
conda activate wgs_env
mkdir -p fastqc_raw   #outputfolder
#run_fastq_on_raw_file
fastqc -t 4 \          
    4-Rest_R1_001.fastq \
    4-Rest_R2_001.fastq \
    -o fastqc_raw

fastqc -t 4 \
    5-Rest_R1_001.fastq \
    5-Rest_R2_001.fastq \
    -o fastqc_raw
    
    fastqc -t 4 \
    5-Susp_R1_001.fastq \
    5-Susp_R2_001.fastq \
    -o fastqc_raw
    
fastqc -t 4 \
    4-Susp_R1_001.fastq \
    4-Susp_R2_001.fastq \
    -o fastqc_raw

#run_fastptrimming_on_raw_file

  mkdir -p fastp_out

R1="4-Rest_R1_001.fastq"
R2="4-Rest_R2_001.fastq"
SAMPLE="4-Rest"

fastp \
  -i "$R1" \
  -I "$R2" \
  -o "fastp_out/${SAMPLE}_R1.trimmed.fastq.gz" \
  -O "fastp_out/${SAMPLE}_R2.trimmed.fastq.gz" \
  --html "fastp_out/${SAMPLE}.html" \
  --json "fastp_out/${SAMPLE}.json" \
  --detect_adapter_for_pe \
  --thread 8

R1="5-Rest_R1_001.fastq"
R2="5-Rest_R2_001.fastq"
SAMPLE="5-Rest"

fastp \
  -i "$R1" \
  -I "$R2" \
  -o "fastp_out/${SAMPLE}_R1.trimmed.fastq.gz" \
  -O "fastp_out/${SAMPLE}_R2.trimmed.fastq.gz" \
  --html "fastp_out/${SAMPLE}.html" \
  --json "fastp_out/${SAMPLE}.json" \
  --detect_adapter_for_pe \
  --thread 8 

R1="4-Susp_R1_001.fastq"
R2="4-Susp_R2_001.fastq"
SAMPLE="4-Susp"

fastp \
  -i "$R1" \
  -I "$R2" \
  -o "fastp_out/${SAMPLE}_R1.trimmed.fastq.gz" \
  -O "fastp_out/${SAMPLE}_R2.trimmed.fastq.gz" \
  --html "fastp_out/${SAMPLE}.html" \
  --json "fastp_out/${SAMPLE}.json" \
  --detect_adapter_for_pe \
  --thread 8


R1="5-Susp_R1_001.fastq"
R2="5-Susp_R2_001.fastq"
SAMPLE="5-Susp"

fastp \
  -i "$R1" \
  -I "$R2" \
  -o "fastp_out/${SAMPLE}_R1.trimmed.fastq.gz" \
  -O "fastp_out/${SAMPLE}_R2.trimmed.fastq.gz" \
  --html "fastp_out/${SAMPLE}.html" \
  --json "fastp_out/${SAMPLE}.json" \
  --detect_adapter_for_pe \
  --thread 8


  #run_fastqc_on_trimmed_file
  mkdir -p fastqc_trimmed
  fastqc fastp_out/*.trimmed.fastq.gz -o fastqc_trimmed -t 8 
  #for_imdividual
  fastqc fastp_out/4-Rest_R1.trimmed.fastq.gz fastp_out/4-Rest_R2.trimmed.fastq.gz \
    -o fastqc_trimmed -t 8

#index_reference_genome
cd /home/g89x126/sc_wgs

# samtools index
samtools faidx reference.fasta

# bwa index
bwa index reference.fasta

#align_trimmed_file_with_reference(align_one.sh)
  #!/bin/bash

conda activate wgs_env

REF="/home/g89x126/sc_wgs/reference.fasta"
TRIMDIR="/home/g89x126/sc_wgs/fastp_out"
OUTDIR="/home/g89x126/sc_wgs/bam"

mkdir -p $OUTDIR

# CHANGE THIS LINE FOR EACH SAMPLE
SAMPLE="4-Rest"

R1="${TRIMDIR}/${SAMPLE}_R1.trimmed.fastq.gz"
R2="${TRIMDIR}/${SAMPLE}_R2.trimmed.fastq.gz"

echo "Aligning sample: $SAMPLE"
echo "Using:"
echo "R1 = $R1"
echo "R2 = $R2"

# ALIGNMENT → BAM SORTING
bwa mem -t 8 $REF $R1 $R2 \
    | samtools view -b -@ 4 - \
    | samtools sort -@ 4 -o ${OUTDIR}/${SAMPLE}.sorted.bam -

# INDEX BAM
samtools index ${OUTDIR}/${SAMPLE}.sorted.bam

echo "DONE: ${OUTDIR}/${SAMPLE}.sorted.bam"


#RUN_VARIANT_CALLING_USING_bcftools
mkdir -p vcf 
bcftools mpileup -Ou -f reference.fasta bam/5-Susp.sorted.bam \
| bcftools call -mv -Oz -o vcf/5-Susp.vcf.gz

#index_vcf_file
bcftools index vcf/5-Susp.vcf.gz



done

