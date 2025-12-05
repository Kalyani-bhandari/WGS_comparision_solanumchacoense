ls ../../
cp ../../snpEff.config .
pwd
nano snpEff.config
snpEff build -gff3 -v solanum_chacoense_m6 -c snpEff.config
ls../solanum_chacoense_m6
cat genes.gff|greo -i disease
