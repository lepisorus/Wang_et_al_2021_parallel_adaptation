#!/bin/bash
INFILE=$1
PREFIX=$(basename $1 | sed 's/\.txt//g')

#file ${PREFIX}.filterRegion.txt contains potential introgression regions from fd analyses and inversion regions (from sliding window PCA analyses)
vcftools --vcf ../../../landrace_palmarChico_filtered.recode.vcf --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --keep ${INFILE} --exclude-bed ${PREFIX}.filterRegion.txt --min-meanDP 10 --maf 0.05 --max-missing-count 0 --recode --out ${PREFIX}

#filter for genotype depth >=10
perl filtersnp.pl ${PREFIX}.recode.vcf > ${PREFIX}.filter.vcf

#filter for SNPs present in each highland vcf files to ensure that the analyses were based on the same set of SNPs
vcftools --vcf ${PREFIX}.filter.vcf --positions finalSNPs.txt --recode --out ${PREFIX}.filter.final 

#calculate allele frequency
vcftools --vcf ${PREFIX}.filter.final.recode.vcf --keep ${PREFIX}2.txt --freq --out ${PREFIX}
vcftools --vcf ${PREFIX}.filter.final.recode.vcf --keep MexLow2.txt --freq --out ${PREFIX}_MexLow
vcftools --vcf ${PREFIX}.filter.final.recode.vcf --keep parv.txt --freq --out ${PREFIX}_parv

#combine together
paste ${PREFIX}.frq ${PREFIX}_MexLow.frq ${PREFIX}_parv.frq | sed 's/\:/\t/g' | sed -e '1d' | cut -f1,2,6,8,14,16,22,24 > ${PREFIX}_MexLow_parv.frq

#calculate Fst
vcftools --vcf ${PREFIX}.filter.final.recode.vcf --weir-fst-pop ${PREFIX}2.txt --weir-fst-pop MexLow2.txt --out ${PREFIX}_MexLow 
vcftools --vcf ${PREFIX}.filter.final.recode.vcf --weir-fst-pop ${PREFIX}2.txt --weir-fst-pop parv.txt --out ${PREFIX}_parv 
vcftools --vcf ${PREFIX}.filter.final.recode.vcf --weir-fst-pop MexLow2.txt --weir-fst-pop parv.txt --out ${PREFIX}_MexLow_parv

#combine Fst result together
paste ${PREFIX}_MexLow.weir.fst ${PREFIX}_parv.weir.fst ${PREFIX}_MexLow_parv.weir.fst | cut -f1,2,3,6,9 > ${PREFIX}.weir.fst


