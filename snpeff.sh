ls ../../
cp ../../snpEff.config .
pwd
nano snpEff.config
snpEff build -gff3 -v solanum_chacoense_m6 -c snpEff.config
ls../solanum_chacoense_m6
cat genes.gff|greo -i disease


#annotate other files with sneff
snpEff -noStats solanum_chacoense_m6 \
   sc_wgs /Susceptible_unique /Susceptible_unique.vcf.gz \
    > EMS.susceptible.annotated.vcf
    
cd /home/g89x126/sc_wgs/snpeff

snpEff -noStats solanum_chacoense_m6 \
    ../Susceptible_unique/Susceptible_unique.vcf.gz \
    > Susceptible.annotated.vcf


snpEff -noStats solanum_chacoense_m6 \
    ../Resistance_unique/Resistance_unique.vcf.gz \
    > Resistance.annotated.vcf

#indexing
bgzip Susceptible.annotated.vcf
bgzip Resistance.annotated.vcf

tabix -p vcf Susceptible.annotated.vcf.gz
tabix -p vcf Resistance.annotated.vcf.gz

#For_Resistance

bcftools view -i 'INFO/ANN ~ "HIGH"' Resistance.annotated.vcf.gz \
    -Oz -o Resistance.HIGH.vcf.gz
    bcftools view -i 'INFO/ANN ~ "HIGH|MODERATE"' Resistance.annotated.vcf.gz \
  -Oz -o Resistance.HIGH.vcf.gz

zcat Resistance.HIGH.vcf.gz | head
zgrep -v "^#" Resistance.HIGH.vcf.gz | head
bcftools view -H Resistance.HIGH.vcf.gz | wc -l

bcftools view -i 'ANN~"MODERATE"' Resistance.annotated.vcf.gz \
  -Oz -o Resistance.MODERATE.vcf.gz
bcftools index Resistance.MODERATE.vcf.gz

bcftools view -H Resistance.MODERATE.vcf.gz | wc -l
zgrep -v "^#" Resistance.MODERATE.vcf.gz | head



bcftools view -i 'INFO/ANN ~ "HIGH"' Susceptible.annotated.vcf.gz \
    -Oz -o Susceptible.HIGH.vcf.gz


    #Extracted_geneID_from_moderate_variants
    zgrep -v "^#" Resistance.MODERATE.vcf.gz \
    | sed 's/.*ANN=//' \
    | tr ',' '\n' \
    | grep "missense_variant" \
    | cut -d '|' -f4 \
    | sort -u \
    > Resistance_MODERATE_genes.txt

    sed -i 's/$/.t1/' Resistance_MODERATE_genes.txt
#alined with functional annotation
sort Resistance_MODERATE_genes.txt > Resistance_MODERATE_genes.sorted.txt
sort functional_annotation.txt > functional_annotation.sorted.txt

sort -k1,1 Resistance_MODERATE_genes.txt > Resistance_MODERATE_genes.sorted.txt
sort -k1,1 functional_annotation.txt > functional_annotation.sorted.txt

join -t $'\t' Resistance_MODERATE_genes.sorted.txt functional_annotation.sorted.txt \
    > Resistance_MODERATE_genes_with_annotation.txt


join Resistance_MODERATE_genes.sorted.txt functional_annotation.sorted.txt \
    > Resistance_MODERATE_genes_with_annotation.txt


#For_Susceptible

bcftools view -i 'INFO/ANN ~ "HIGH"' Resistance.annotated.vcf.gz \
    -Oz -o Resistance.HIGH.vcf.gz
    
    bcftools view -i 'INFO/ANN ~ "HIGH|MODERATE"' Susceptible.annotated.vcf.gz \
  -Oz -o Susceptible.HIGH.vcf.gz
  bcftools index Susceptible.HIGH.vcf.gz

zcat Susceptible.HIGH.vcf.gz | head
zgrep -v "^#" Susceptible.HIGH.vcf.gz | head
bcftools view -H Susceptible.HIGH.vcf.gz | wc -l

bcftools view -i 'ANN~"MODERATE"' Susceptible.annotated.vcf.gz \
  -Oz -o Susceptible.MODERATE.vcf.gz
bcftools index Susceptible.MODERATE.vcf.gz

bcftools view -H Susceptible.MODERATE.vcf.gz | wc -l
zgrep -v "^#" Susceptible.MODERATE.vcf.gz | head



bcftools view -i 'INFO/ANN ~ "HIGH"' Susceptible.annotated.vcf.gz \
    -Oz -o Susceptible.HIGH.vcf.gz


    #Extracted_geneID_from_moderate_variants
    zgrep -v "^#" Susceptible.MODERATE.vcf.gz \
    | sed 's/.*ANN=//' \
    | tr ',' '\n' \
    | grep "missense_variant" \
    | cut -d '|' -f4 \
    | sort -u \
    > Susceptible_MODERATE_genes.txt

    sed -i 's/$/.t1/' Susceptible_MODERATE_genes.txt
#alined with functional annotation

sort -k1,1 Susceptible_MODERATE_genes.txt > Susceptible_MODERATE_genes.sorted.txt
sort -k1,1 functional_annotation.txt > functional_annotation.sorted.txt

join -t $'\t' Susceptible_MODERATE_genes.sorted.txt functional_annotation.sorted.txt \
    > Susceptible_MODERATE_genes_with_annotation.txt





