#!/bin/bash
#SBATCH --job-name=align
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=nextgen
#SBATCH --output=align_%A_%a.out
#SBATCH --error=align_%A_%a.err

module load bwa
module load samtools

RAW_DIR="/home/g89x126/sc_wgs/fastp_out"
REF="/home/g89x126/sc_wgs/reference.fasta"
OUTDIR="/home/g89x126/sc_wgs/bam"

mkdir -p $OUTDIR logs

# ---------------------------------------------------------
# 3. Create sample list automatically
# ---------------------------------------------------------
SAMPLES=( $(ls $RAW_DIR/*_R1_trimmed.fastq.gz | sed 's/_R1_trimmed.fastq.gz//' | xargs -n1 basename) )

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

R1="$RAW_DIR/${SAMPLE}_R1_trimmed.fastq.gz"
R2="$RAW_DIR/${SAMPLE}_R2_trimmed.fastq.gz"

echo "Processing sample: $SAMPLE"
echo "R1 = $R1"
echo "R2 = $R2"

# ---------------------------------------------------------
# 4. Index reference (done only once)
# ---------------------------------------------------------
if [ ! -f "${REF}.bwt" ]; then
    echo "Indexing reference genome with BWA..."
    bwa index $REF
fi

# ---------------------------------------------------------
# 5. Alignment → Sorting → BAM index
# ---------------------------------------------------------
bwa mem -t 8 $REF $R1 $R2 \
    | samtools view -bS -@ 4 - \
    | samtools sort -@ 4 -o $OUTDIR/${SAMPLE}.sorted.bam -

samtools index $OUTDIR/${SAMPLE}.sorted.bam

echo "Finished sample: $SAMPLE"