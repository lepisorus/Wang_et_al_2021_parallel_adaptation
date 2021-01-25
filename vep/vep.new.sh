#!/bin/bash

PREFIX=$1 # the pair of highland populations


vcftools --vcf ../MexHigh.filter.final.recode.vcf --positions ../${PREFIX}.common.outlierPBE.txt --recode --out ${PREFIX}.common.outlierPBE
cat ${PREFIX}.common.outlierPBE.recode.vcf |  grep -v "^#" | awk -v OFS="\t" '{print $1, $2, $2, $4"/"$5, "+"}' > ${PREFIX}.forVep
rm ${PREFIX}.common.outlierPBE.recode.vcf
/work/LAS/mhufford-lab/lwang/bin/ensembl-vep/vep --dir /work/LAS/mhufford-lab/lwang/bin/ensembl-vep/zea_mays_cache --species zea_mays --cache_version 80 -i ${PREFIX}.forVep --output_file ${PREFIX}.Vepout --stats_file ${PREFIX}.forVep.stats --stats_text --offline --hgvs --force_overwrite --fasta /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.22.dna.genome.fa

