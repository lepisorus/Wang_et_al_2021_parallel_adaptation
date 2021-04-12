#!/bin/bash
INFILE=$1
PREFIX=$(basename $1 | sed 's/\.filter\.final\.recode\.vcf//g')

awk -v OFS="\t" '{print $1"_"$2}' ${INFILE} > ${PREFIX}.ID
paste <(cut -f1,2 ${INFILE})  ${PREFIX}.ID <(cut -f4- ${INFILE}) | sed 's/ /\t/g' > ${PREFIX}.4thinning.vcf
 
plink --vcf ${PREFIX}.4thinning.vcf --make-bed --out ${PREFIX}.test

plink --bfile ${PREFIX}.test --indep-pairwise 10 1 0.2 --out ${PREFIX}

sed -i 's/\_/\t/g' ${PREFIX}.prune.in

vcftools --vcf ${PREFIX}.4thinning.vcf --positions  ${PREFIX}.prune.in --recode --out ${PREFIX}.LDprune
 

#calculate Fst
vcftools --vcf ${PREFIX}.LDprune.recode.vcf --weir-fst-pop ${PREFIX}2.txt --weir-fst-pop MexLow2.txt --out ${PREFIX}_MexLow.prune 
vcftools --vcf ${PREFIX}.LDprune.recode.vcf --weir-fst-pop ${PREFIX}2.txt --weir-fst-pop parv.txt --out ${PREFIX}_parv.prune 
vcftools --vcf ${PREFIX}.LDprune.recode.vcf --weir-fst-pop MexLow2.txt --weir-fst-pop parv.txt --out ${PREFIX}_MexLow_parv.prune

paste ${PREFIX}_MexLow.prune.weir.fst ${PREFIX}_parv.prune.weir.fst ${PREFIX}_MexLow_parv.prune.weir.fst | cut -f1,2,3,6,9 > ${PREFIX}.weir.fst2

Rscript calcPBE.R -i ${PREFIX}.weir.fst2 -o ${PREFIX}.allPBE.prune.txt -s ${PREFIX}.outlierPBE.prune.txt -n ${PREFIX}.neutralPBE.prune.txt