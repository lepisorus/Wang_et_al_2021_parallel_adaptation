bedtools intersect -a <(grep -v "^CHROM" MexHigh.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > MexHighOutlier.withinGene.txt
bedtools intersect -a <(grep -v "^CHROM" MexHigh.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > MexHighOutlier.emptySNP.txt
bedtools closest -a MexHighOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > MexHighOutlier.closestGene.txt
cat <(cut -f7 MexHighOutlier.withinGene.txt) <(cut -f7 MexHighOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > MexHighOutlier.fullGene.txt

bedtools intersect -a <(grep -v "^CHROM" GuaHigh.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > GuaHighOutlier.withinGene.txt
bedtools intersect -a <(grep -v "^CHROM" GuaHigh.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > GuaHighOutlier.emptySNP.txt
bedtools closest -a GuaHighOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > GuaHighOutlier.closestGene.txt
cat <(cut -f7 GuaHighOutlier.withinGene.txt) <(cut -f7 GuaHighOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > GuaHighOutlier.fullGene.txt

bedtools intersect -a <(grep -v "^CHROM" Andes.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > AndesOutlier.withinGene.txt
bedtools intersect -a <(grep -v "^CHROM" Andes.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > AndesOutlier.emptySNP.txt
bedtools closest -a AndesOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > AndesOutlier.closestGene.txt
cat <(cut -f7 AndesOutlier.withinGene.txt) <(cut -f7 AndesOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > AndesOutlier.fullGene.txt


bedtools intersect -a <(grep -v "^CHROM" SW_US.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -wa -wb > SW_USOutlier.withinGene.txt
bedtools intersect -a <(grep -v "^CHROM" SW_US.outlierPBE.prune.txt | awk -v OFS="\t" '{print $1, $2, $2}' | sort -n -k1,1 -n -k2,2 ) -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -v > SW_USOutlier.emptySNP.txt
bedtools closest -a SW_USOutlier.emptySNP.txt -b /work/LAS/mhufford-lab/lwang/AGPv3/zea_mays.protein_coding.bed -D ref | awk -v OFS="\t" '$8<=10000 && $8 >= -10000' > SW_USOutlier.closestGene.txt
cat <(cut -f7 SW_USOutlier.withinGene.txt) <(cut -f7 SW_USOutlier.closestGene.txt) | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > SW_USOutlier.fullGene.txt


