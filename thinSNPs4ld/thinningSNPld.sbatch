#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=96:00:00
#SBATCH --job-name=R
#SBATCH --output=R_%j.out 
#SBATCH --error=R_%j.err
#SBATCH --mail-user=lilepisorus@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module load plink
module load vcftools
module load r
module load r-data-table


./thinningSNPld2.sh MexHigh.filter.final.recode.vcf
./thinningSNPld2.sh GuaHigh.filter.final.recode.vcf
./thinningSNPld2.sh SW_US.filter.final.recode.vcf

vcftools --vcf Andes.LDprune.recode.vcf --weir-fst-pop Andes2.txt --weir-fst-pop SA_Low2.txt --out Andes_SA_Low.prune 
vcftools --vcf Andes.LDprune.recode.vcf --weir-fst-pop Andes2.txt --weir-fst-pop parv.txt --out Andes_parv.prune  
vcftools --vcf Andes.LDprune.recode.vcf --weir-fst-pop SA_Low2.txt --weir-fst-pop parv.txt --out Andes_SA_Low_parv.prune 

paste Andes_SA_Low.prune.weir.fst Andes_parv.prune.weir.fst Andes_SA_Low_parv.prune.weir.fst | cut -f1,2,3,6,9 > Andes.weir.fst2
Rscript calcPBE.R -i Andes.weir.fst2 -o Andes.allPBE.prune.txt -s Andes.outlierPBE.prune.txt -n Andes.neutralPBE.prune.txt

