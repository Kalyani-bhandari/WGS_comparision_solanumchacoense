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
  
done

