#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > MexHigh.Outlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../SW_US.outlierPBE.txt | sort -n -k1,1 -n -k2,2) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > SW_US.Outlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.outlierPBE.txt | sort -n -k1,1 -n -k2,2) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > Andes.Outlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../GuaHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > GuaHigh.Outlier.annotation.txt

##bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > MexHighOutlier.withinGene.txt
##bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > MexHighOutlier.emptySNP.txt
#bedtools closest -a MexHighOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D b > MexHighOutlier.closestGene.txt
#cat <(cut -f7 MexHighOutlier.withinGene.txt) <(cut -f7 MexHighOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHighOutlier.fullGene.txt

#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../SW_US.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > SW_USOutlier.withinGene.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../SW_US.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > SW_USOutlier.emptySNP.txt
#bedtools closest -a SW_USOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > SW_USOutlier.closestGene.txt
#cat <(cut -f7 SW_USOutlier.withinGene.txt) <(cut -f7 SW_USOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > SW_USOutlier.fullGene.txt

#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > Andes.Outlier.withinGene.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > Andes.Outlier.emptySNP.txt
#bedtools closest -a Andes.Outlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > Andes.Outlier.closestGene.txt
#cat <(cut -f7 Andes.Outlier.withinGene.txt) <(cut -f7 Andes.Outlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > Andes.Outlier.fullGene.txt

#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../GuaHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > GuaHigh.Outlier.withinGene.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../GuaHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > GuaHigh.Outlier.emptySNP.txt
#bedtools closest -a GuaHigh.Outlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > GuaHigh.Outlier.closestGene.txt
#cat <(cut -f7 GuaHigh.Outlier.withinGene.txt) <(cut -f7 GuaHigh.Outlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > GuaHigh.Outlier.fullGene.txt

bedtools closest -a MexHighOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > MexHighOutlier.10kbGene.txt
cat <(cut -f7 MexHighOutlier.withinGene.txt) <(cut -f7 MexHighOutlier.10kbGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHighOutlier.fullGene.txt

bedtools closest -a GuaHigh.Outlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > GuaHigh.Outlier.10kbGene.txt
cat <(cut -f7 GuaHigh.Outlier.withinGene.txt) <(cut -f7 GuaHigh.Outlier.10kbGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > GuaHigh.Outlier.fullGene.txt

bedtools closest -a SW_USOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > SW_USOutlier.10kbGene.txt
cat <(cut -f7 SW_USOutlier.withinGene.txt) <(cut -f7 SW_USOutlier.10kbGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > SW_USOutlier.fullGene.txt

bedtools closest -a Andes.Outlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > Andes.Outlier.10kbGene.txt
cat <(cut -f7 Andes.Outlier.withinGene.txt) <(cut -f7 Andes.Outlier.10kbGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > Andes.Outlier.fullGene.txt

bedtools intersect -a ../finalSNPs.bed -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > finalSNPs.withinGene.txt
bedtools intersect -a ../finalSNPs.bed -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > finalSNPs.emptySNP.txt
bedtools closest -a finalSNPs.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > finalSNPs.10kbGene.txt
cat <(cut -f7 finalSNPs.withinGene.txt) <(cut -f7 finalSNPs.10kbGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > finalSNPs.fullGene.txt
