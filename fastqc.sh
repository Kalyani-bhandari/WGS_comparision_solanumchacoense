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
samtools faidx reference.fa

# bwa index
bwa index reference.fa

#align_trimmed_file_with_reference(align_one.sh)
  #!/bin/bash

conda activate wgs_env

REF="/home/g89x126/sc_wgs/reference.fa"
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

#filter_vcf_file
mkdir -p vcf_filtered


bcftools filter \
   -e 'QUAL<30 || DP<8 || MQ<30 || MQ0F>0.05 || VDB<0.0001 || SGB<-0.6' \
   vcf/4-Rest.vcf.gz -Oz -o vcf_filtered/4-Rest.filtered.vcf.gz

   bcftools index vcf_filtered/4-Rest.filtered.vcf.gz


bcftools filter \
   -e  'QUAL<30 || DP<8 || MQ<30 || MQ0F>0.05 || VDB<0.0001 || SGB<-0.6' \
   vcf/4-Susp.vcf.gz -Oz -o vcf_filtered/4-Susp.filtered.vcf.gz

bcftools index vcf_filtered/4-Susp.filtered.vcf.gz

bcftools filter \
   -e  'QUAL<30 || DP<8 || MQ<30 || MQ0F>0.05 || VDB<0.0001 || SGB<-0.6' \
   vcf/5-Susp.vcf.gz -Oz -o vcf_filtered/5-Susp.filtered.vcf.gz

bcftools index vcf_filtered/5-Susp.filtered.vcf.gz

bcftools filter \
   -e  'QUAL<30 || DP<8 || MQ<30 || MQ0F>0.05 || VDB<0.0001 || SGB<-0.6' \
   vcf/5-Rest.vcf.gz -Oz -o vcf_filtered/5-Rest.filtered.vcf.gz

bcftools index vcf_filtered/5-Rest.filtered.vcf.gz

#Extract_EMS_Mutation_only
mkdir -p EMS_vcfs

for VCF in vcf_filtered/*.filtered.vcf.gz; do
    BASE=$(basename "$VCF" .filtered.vcf.gz)
    
    bcftools view \
      -i 'TYPE="snp" && ((REF="G" && ALT="A") || (REF="C" && ALT="T"))' \
      "$VCF" \
      -Oz -o EMS_vcfs/${BASE}.EMS.vcf.gz
    
    bcftools index EMS_vcfs/${BASE}.EMS.vcf.gz
done
#countit
bcftools stats EMS_vcfs/4-Rest.EMS.vcf.gz | grep "number of records"


#FILTER_ONLY_REST_OR_ONLY_SUSP

bcftools isec \
    -p compare4 \
    EMS_vcfs/4-Rest.EMS.vcf.gz \
    EMS_vcfs/4-Susp.EMS.vcf.gz


bgzip compare4/0000.vcf
bcftools index compare4/0000.vcf.gz

bgzip compare4/0001.vcf
bcftools index compare4/0001.vcf.gz

bcftools isec \
    -p compare5 \
    EMS_vcfs/5-Rest.EMS.vcf.gz \
    EMS_vcfs/5-Susp.EMS.vcf.gz

bgzip compare5/0000.vcf
bcftools index compare5/0000.vcf.gz

bgzip compare5/0001.vcf
bcftools index compare5/0001.vcf.gz

#count
bcftools view -H compare4/4-Rest_only.vcf.gz | wc -l    # unique to 4-Rest
bcftools view -H compare4/4-Susp_only.vcf.gz | wc -l    # unique to 4-Susp
bcftools view -H compare4/0002.vcf.gz| wc -l    # shared

bcftools view -H compare4/4-Rest_only.vcf.gz | head
bcftools view -H compare4/4-Susp_only.vcf.gz | head


bcftools view -H compare4/4-Rest_only.vcf.gz | wc -l
bcftools view -H compare4/4-Susp_only.vcf.gz | wc -l
bcftools view -H compare4/0002.vcf.gz | wc -l


#common snps in both resistance lines 4 and that is not in susp

# 1. Shared resistant SNPs
bcftools isec -n=2 EMS_vcfs/4-Rest.EMS.vcf.gz EMS_vcfs/5-Rest.EMS.vcf.gz -p Rest_shared
bgzip Rest_shared/0000.vcf
bcftools index Rest_shared/0000.vcf.gz

# 2. Shared susceptible SNPs
bcftools isec -n=2 EMS_vcfs/4-Susp.EMS.vcf.gz EMS_vcfs/5-Susp.EMS.vcf.gz -p Susp_shared

bgzip Susp_shared/0000.vcf
bcftools index Susp_shared/0000.vcf.gz

# 3. Resistance-specific SNPs
bcftools isec Rest_shared/0000.vcf.gz Susp_shared/0000.vcf.gz -p Resistance_unique

bgzip Resistance_unique/0000.vcf
bcftools index Resistance_unique/0000.vcf.gz

# 4. Count final candidates
bcftools view -H Resistance_unique/0000.vcf.gz | wc -l
bcftools view -H Resistance_unique/Resistance_unique.vcf.gz | head

#5. Susceptible_unique_SNPs
bcftools isec \
    Susp_shared/0000.vcf.gz \
    Rest_shared/0000.vcf.gz \
    -p Susceptible_unique

bgzip Susceptible_unique/0000.vcf
bcftools index Susceptible_unique/0000.vcf.gz


# 6. Count final candidates
bcftools view -H Susceptible_unique/0000.vcf.gz | wc -l















#annotation
#snepff_database
mkdir -p snpeff/data/solanum_chacoense_m6

cp reference.fa snpeff/data/solanum_chacoense_m6/sequences.fa
cp m6.hc.gene_models.all_chrs_corrected.gff3 snpeff/data/solanum_chacoense_m6/genes.gff

find snpeff -name snpEff.config
nano snpeff/snpEff.config
solanum_chacoense_m6.genome : solanum_chacoense_m6

cd snpeff
java -Xmx8g -jar snpEff.jar build -gff3 -v solanum_chacoense_m6

cd snpeff

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

mkdir -p $CONDA_PREFIX/share/snpeff/data/solanum_chacoense_m6

conda activate wgs_env
conda install -c bioconda snpeff

bcftools csq \
   --force \
   -f reference.fa \
   -g m6.hc.gene_models.all_chrs_corrected.gff3 \
   EMS_vcfs/4-Rest.EMS.vcf.gz \
   -Oz -o EMS_vcfs/4-Rest.EMS.annotated.vcf.gz

   #another_file_for_annotation
   bcftools csq --force \
    -f reference.fa \
    -g m6.hc.gene_models.chr_corrected.gff3 \
    EMS_vcfs/4-Rest.EMS.vcf.gz \
    -Oz -o EMS_vcfs/4-Rest.EMS.annotated.v2.vcf.gz

   
bcftools +split-vep EMS_vcfs/4-Rest.EMS.annotated.vcf.gz \
    -d both \
    -Oz -o EMS_vcfs/4-Rest.EMS.annotated.split.vcf.gz
bcftools +split-vep -l -a BCSQ EMS_vcfs/4-Rest.EMS.annotated.vcf.gz

bcftools +split-vep \
    -a BCSQ \
    -c Consequence,gene,transcript,amino_acid_change,dna_change \
    -s worst \
    -p BCSQ \
    EMS_vcfs/4-Rest.EMS.annotated.vcf.gz \
    -Oz -o EMS_vcfs/4-Rest.EMS.split.vcf.gz

bcftools index -f EMS_vcfs/4-Rest.EMS.split.vcf.gz




#Compare_Resistant_vs_Susceptible_VCFs

#Resistant_shared_mutations:
mkdir -p compare/resistant

bcftools isec \
   vcf_filtered/4-Rest.filtered.vcf.gz \
   vcf_filtered/5-Rest.filtered.vcf.gz \
   -p compare/resistant
   bgzip compare/resistant/0002.vcf
tabix -p vcf compare/resistant/0002.vcf.gz


#Susceptible_shared_mutations:
mkdir -p compare/susceptible

bcftools isec \
   vcf_filtered/4-Susp.filtered.vcf.gz \
   vcf_filtered/5-Susp.filtered.vcf.gz \
   -p compare/susceptible
   bgzip compare/susceptible/0002.vcf
tabix -p vcf compare/susceptible/0002.vcf.gz

#Mutation_only_in_resistant
mkdir -p compare/final

bcftools isec \
   compare/resistant/0002.vcf.gz \
   compare/susceptible/0002.vcf.gz \
   -p compare/final
#indexingandcompress
bgzip compare/final/0000.vcf
tabix -p vcf compare/final/0000.vcf.gz

bgzip compare/final/0001.vcf
tabix -p vcf compare/final/0001.vcf.gz

bgzip compare/final/0002.vcf
tabix -p vcf compare/final/0002.vcf.gz

#Filter EMS-type SNPs from final/0000.vcf.gz
mkdir -p compare/ems_filtered

bcftools view \
   -i '(REF="G" && ALT="A") || (REF="C" && ALT="T")' \
   compare/final/0000.vcf.gz \
   -Oz -o compare/ems_filtered/0000.EMS.vcf.gz
#index_and_compress
tabix -p vcf compare/ems_filtered/0000.EMS.vcf.gz

#Count_SNPs_and_Indels
bcftools view -H compare/final/0000.vcf.gz | wc -l
bcftools view -H compare/ems_filtered/0000.EMS.vcf.gz | wc -l
bcftools view -i 'TYPE="indel"' compare/final/0000.vcf.gz | wc -l

#Extract_high-impact_mutations
bcftools filter -i 'MQ>=40 && DP>=10 && TYPE="snp"' \
   compare/ems_filtered/0000.EMS.vcf.gz \
| bcftools csq --force \
   -f reference.fa \
   -g reference.gff3 \
   -Oz -o EMS.csq.vcf.gz
grep -v "#" reference.gff3 | cut -f1 | head
grep -v "#" m6.hc.gene_models.chr_corrected.gff3 | cut -f1 | head

bcftools csq --force \
   -f reference.fa \
   -g m6.hc.gene_models.chr_corrected.gff3 \
   -Oz -o EMS.annotated.vcf.gz \
   compare/ems_filtered/0000.EMS.vcf.gz




done

