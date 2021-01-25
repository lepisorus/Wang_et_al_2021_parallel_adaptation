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
#module load java
module load vcftools
module load jdk/10.0.2_13-fr57jru

vcftools --vcf ../../../../landrace_palmarChico_filtered.recode.vcf --positions ../finalSNPs.txt --recode --out landrace_palmarChico_PBEfinalSNP_v3
java -Xmx120g -jar ./bin/snpEff/snpEff.jar -c ./bin/snpEff/snpEff.config Zea_mays_AGPv3 landrace_palmarChico_PBEfinalSNP_v3.recode.vcf > landrace_palmarChico_PBEfinalSNP_v3.SNPeff.vcf

cat <(awk '{print $1"_"$2}' ../MexHigh.outlierPBE) <(awk '{print $1"_"$2}' ../GuaHigh.outlierPBE) <(awk '{print $1"_"$2}' ../SW_US.outlierPBE) <(awk '{print $1"_"$2}' ../Andes.outlierPBE) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > allPop.outlierSites
vcftools --vcf /work/LAS/mhufford-lab/lwang/landrace_palmarChico_filtered.recode.vcf --positions allPop.outlierSites --recode --out allPop.outlierSites
java -Xmx120g -jar ./bin/snpEff/snpEff.jar -c ./bin/snpEff/snpEff.config Zea_mays_AGPv3 allPop.outlierSites.recode.vcf > allPop.outlierSites.SNPeff.vcf

#ssh condo "qstat -f ${SLURM_JOB_ID} |head"

less allPop.outlierSites.SNPeff.vcf | grep -v "^#" | sed 's/|/\t/g' | awk '{print $9}' > allPop.outlierSites.SNPeff.txt

