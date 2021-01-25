#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' MexHigh.GuaHigh.SW_US.commonOutlier.txt) -b ../Zea_mays.AGPv3.21.bed -wo > MexHigh.GuaHigh.SW_US.commonOutlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' MexHigh.GuaHigh.SW_US.Andes.commonOutlier.txt) -b ../Zea_mays.AGPv3.21.bed -wo > MexHigh.GuaHigh.SW_US.Andes.commonOutlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > MexHigh.GuaHigh.commonOutlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.SW_US.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > MexHigh.SW_US.commonOutlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.Andes.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > MexHigh.Andes.commonOutlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > Andes.GuaHigh.commonOutlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../SW_US.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > SW_US.GuaHigh.commonOutlier.annotation.txt
#bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.SW_US.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo > Andes.SW_US.commonOutlier.annotation.txt

bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > MexHigh.GuaHigh.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > MexHigh.GuaHigh.commonOutlier.emptySNP.txt
bedtools closest -a MexHigh.GuaHigh.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > MexHigh.GuaHigh.commonOutlier.closestGene.txt
cat <(cut -f7 MexHigh.GuaHigh.commonOutlier.withinGene.txt) <(cut -f7 MexHigh.GuaHigh.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHigh.GuaHigh.commonOutlier.fullGene.txt

bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.SW_US.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > MexHigh.SW_US.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.SW_US.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > MexHigh.SW_US.commonOutlier.emptySNP.txt
bedtools closest -a MexHigh.SW_US.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > MexHigh.SW_US.commonOutlier.closestGene.txt
cat <(cut -f7 MexHigh.SW_US.commonOutlier.withinGene.txt) <(cut -f7 MexHigh.SW_US.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHigh.SW_US.commonOutlier.fullGene.txt

bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.Andes.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > MexHigh.Andes.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../MexHigh.Andes.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > MexHigh.Andes.commonOutlier.emptySNP.txt
bedtools closest -a MexHigh.Andes.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > MexHigh.Andes.commonOutlier.closestGene.txt
cat <(cut -f7 MexHigh.Andes.commonOutlier.withinGene.txt) <(cut -f7 MexHigh.Andes.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHigh.Andes.commonOutlier.fullGene.txt

bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > Andes.GuaHigh.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > Andes.GuaHigh.commonOutlier.emptySNP.txt
bedtools closest -a Andes.GuaHigh.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > Andes.GuaHigh.commonOutlier.closestGene.txt
cat <(cut -f7 Andes.GuaHigh.commonOutlier.withinGene.txt) <(cut -f7 Andes.GuaHigh.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > Andes.GuaHigh.commonOutlier.fullGene.txt

bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../SW_US.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > SW_US.GuaHigh.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../SW_US.GuaHigh.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > SW_US.GuaHigh.commonOutlier.emptySNP.txt
bedtools closest -a SW_US.GuaHigh.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > SW_US.GuaHigh.commonOutlier.closestGene.txt
cat <(cut -f7 SW_US.GuaHigh.commonOutlier.withinGene.txt) <(cut -f7 SW_US.GuaHigh.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > SW_US.GuaHigh.commonOutlier.fullGene.txt

bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.SW_US.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > Andes.SW_US.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' ../Andes.SW_US.common.outlierPBE.txt) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > Andes.SW_US.commonOutlier.emptySNP.txt
bedtools closest -a Andes.SW_US.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > Andes.SW_US.commonOutlier.closestGene.txt
cat <(cut -f7 Andes.SW_US.commonOutlier.withinGene.txt) <(cut -f7 Andes.SW_US.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > Andes.SW_US.commonOutlier.fullGene.txt

cat <(sed '1d' ../MexHigh.outlierPBE.txt) <(sed '1d' ../GuaHigh.outlierPBE.txt) <(sed '1d' ../SW_US.outlierPBE.txt) | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==3{print $2, $3}' > MexHigh.GuaHigh.SW_US.commonOutlier
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' MexHigh.GuaHigh.SW_US.commonOutlier) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > MexHigh.GuaHigh.SW_US.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' MexHigh.GuaHigh.SW_US.commonOutlier) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > MexHigh.GuaHigh.SW_US.commonOutlier.emptySNP.txt
bedtools closest -a MexHigh.GuaHigh.SW_US.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > MexHigh.GuaHigh.SW_US.commonOutlier.closestGene.txt
cat <(cut -f7 MexHigh.GuaHigh.SW_US.commonOutlier.withinGene.txt) <(cut -f7 MexHigh.GuaHigh.SW_US.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHigh.GuaHigh.SW_US.commonOutlier.fullGene.txt

cat <(sed '1d' ../MexHigh.outlierPBE.txt) <(sed '1d' ../GuaHigh.outlierPBE.txt) <(sed '1d' ../SW_US.outlierPBE.txt)  <(sed '1d' ../Andes.outlierPBE.txt) | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==4{print $2, $3}' > MexHigh.GuaHigh.SW_US.Andes.commonOutlier
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' MexHigh.GuaHigh.SW_US.Andes.commonOutlier) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > MexHigh.GuaHigh.SW_US.Andes.commonOutlier.withinGene.txt
bedtools intersect -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2}' MexHigh.GuaHigh.SW_US.Andes.commonOutlier) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > MexHigh.GuaHigh.SW_US.Andes.commonOutlier.emptySNP.txt
bedtools closest -a MexHigh.GuaHigh.SW_US.Andes.commonOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed > MexHigh.GuaHigh.SW_US.Andes.commonOutlier.closestGene.txt
cat <(cut -f7 MexHigh.GuaHigh.SW_US.Andes.commonOutlier.withinGene.txt) <(cut -f7 MexHigh.GuaHigh.SW_US.Andes.commonOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHigh.GuaHigh.SW_US.Andes.commonOutlier.fullGene.txt
