# To Do Steps
Some basic knowledge: [https://github.com/dpaudel/NepalBioinformatics]

1. Run FASTQC on raw files 
2. Run fastp on raw files
3. Run FASTQC on trimmed files
4. Index reference.fasta using ```samtools faidx```
5. Run ```bwa``` for alingment to reference.fasta > output .sam file
6. Sort ```.sam``` file using ```samtools sort```
7. Convert sorted .sam file to bam file
8. Index bam file ```samtools index```
9. Install IGV Integrated Genome Viewer [https://igv.org/].
10. Download ```.bam``` file and ```.bai``` index and ```reference.fasta``` and ```reference.fasta.fai```
11. View the alignment in IGV viewer
