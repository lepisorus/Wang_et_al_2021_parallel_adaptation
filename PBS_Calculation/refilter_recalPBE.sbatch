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
module load vcftools
module load bedtools2
module load R
module load perl

##run the first two commands in the vcftools.sh, filter the SNPs for the described criterion in each population
#vcftools --vcf ../../../landrace_palmarChico_filtered.recode.vcf --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --keep ${INFILE} --exclude-bed ${PREFIX}.filterRegion.txt --min-meanDP 10 --maf 0.05 --max-missing-count 0 --recode --out ${PREFIX}
#perl filtersnp.pl ${PREFIX}.recode.vcf > ${PREFIX}.filter.vcf
./vcftools.sh GuaHigh.txt
./vcftools.sh MexHigh.txt
./vcftools.sh SW_US.txt
./vcftools_Andes.sh Andes.txt

#find the common SNP set for all four highland-lowland-parv vcfs
cat <(less MexHigh.filter.vcf | grep -v "^" | awk -v OFS="" '{print $1, "_", $2}') <(less SW_US.filter.vcf | grep -v "^" | awk -v OFS="" '{print $1, "_", $2}') <(less GuaHigh.filter.vcf | grep -v "^" | awk -v OFS="" '{print $1, "_", $2}') <(less Andes.filter.vcf | grep -v "^" | awk -v OFS="" '{print $1, "_", $2}') | sort -k1,1 | uniq -c | sed 's/^[\t ]*//g' | awk -v FS=" " '$1==4 {print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > finalSNPs.txt

#run the remaining command in vcftools.sh except the first two commands
./vcftools.sh GuaHigh.txt
./vcftools.sh MexHigh.txt
./vcftools.sh SW_US.txt
./vcftools_Andes.sh Andes.txt

#calculate PBE values for each site
Rscript calcPBE.R -i MexHigh.weir.fst -o MexHigh.allPBE.txt -s MexHigh.outlierPBE.txt -n MexHigh.neutralPBE.txt
Rscript calcPBE.R -i GuaHigh.weir.fst -o GuaHigh.allPBE.txt -s GuaHigh.outlierPBE.txt -n GuaHigh.neutralPBE.txt
Rscript calcPBE.R -i SW_US.weir.fst -o SW_US.allPBE.txt -s SW_US.outlierPBE.txt -n SW_US.neutralPBE.txt
Rscript calcPBE.R -i Andes.weir.fst -o Andes.allPBE.txt -s Andes.outlierPBE.txt -n Andes.neutralPBE.txt

#determine if the change of SNP frequency is co-directional, anti-directional, or no-directional for each pairwise highland comparison 
R CMD BATCH determineSNPtype.R

#find out the common outlier and common neutral sites for each pairwise highland comparison
./commonOutlierNeutral.sh
#summarize SNP frequency change type in outlier SNPs and neutral SNPs
R CMD BATCH summarySNPtype.R
#bar plot for co-/anti- direction SNPs in all six pairs of highland populations and for neutral and outlier SNPs
R CMD BATCH co_anti_dir_summary_barPlot.R 




