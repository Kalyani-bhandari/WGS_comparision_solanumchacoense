conda activate wgs_env
mkdir -p fastqc_raw   #outputfolder
fastqc -t 4 \          #run fastqc
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

    mkdir -p fastp_out  #fastp trimming
    for R1 in *_R1_001.fastq; do
    R2="${R1/_R1_001.fastq/_R2_001.fastq}"
    SAMPLE=$(basename "$R1" _R1_001.fastq)

    echo "Running fastp on sample: $SAMPLE"

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "fastp_out/${SAMPLE}_R1.trimmed.fastq.gz" \
        -O "fastp_out/${SAMPLE}_R2.trimmed.fastq.gz" \
        --html "fastp_out/${SAMPLE}.html" \
        --json "fastp_out/${SAMPLE}.json" \
        --detect_adapter_for_pe \
        --thread 8
done

