cat <(sed '1d' SW_US.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') \
<(sed '1d' MexHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') \
<(sed '1d' GuaHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') \
<(sed '1d' Andes.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') \
| sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==4{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/fourHighOutlierSNP.txt

cat <(sed '1d' SW_US.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') \
<(sed '1d' MexHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') \
<(sed '1d' GuaHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') \
| sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==3{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/threeHighMesoAmeOutlierSNP.txt

cat <(sed '1d' MexHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' GuaHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/MH.GH.commonOutlierSNP.txt
cat <(sed '1d' MexHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' SW_US.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/MH.US.commonOutlierSNP.txt
cat <(sed '1d' MexHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' Andes.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/MH.AN.commonOutlierSNP.txt
cat <(sed '1d' SW_US.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' Andes.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/US.AN.commonOutlierSNP.txt
cat <(sed '1d' GuaHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' Andes.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/GH.AN.commonOutlierSNP.txt
cat <(sed '1d' GuaHigh.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' SW_US.outlierPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./outlierSNPs/GH.US.commonOutlierSNP.txt

cat <(sed '1d' MexHigh.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' GuaHigh.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./MH.GH.commonNeutralSNP.txt
cat <(sed '1d' MexHigh.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' SW_US.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./MH.US.commonNeutralSNP.txt
cat <(sed '1d' MexHigh.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' Andes.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./MH.AN.commonNeutralSNP.txt
cat <(sed '1d' GuaHigh.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' Andes.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./GH.AN.commonNeutralSNP.txt
cat <(sed '1d' SW_US.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' Andes.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./US.AN.commonNeutralSNP.txt
cat <(sed '1d' GuaHigh.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') <(sed '1d' SW_US.neutralPBE.txt | awk -v OFS="\t" '{print $1"_"$2}') | sort -k1,1 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2}' | sed 's/\_/\t/g' | sort -n -k1,1 -n -k2,2 > ./GH.US.commonNeutralSNP.txt


vcftools --vcf landrace_palmarChico_PBEfinalSNP_v3.recode.vcf --positions MH.GH.commonOutlierSNP.txt --thin 2000 --recode --out PBEfinalSNP.MH.GH.Outlier.thin2kb
less PBEfinalSNP.MH.GH.Outlier.thin2kb.recode.vcf | grep -v '^#' | awk -v OFS="\t" '{start=$2-10000; end=$2+10000; print $1, start, end}' > MH.GH.OutlierSNP.focalRegions.txt

vcftools --vcf landrace_palmarChico_PBEfinalSNP_v3.recode.vcf --positions MH.US.commonOutlierSNP.txt --thin 2000 --recode --out PBEfinalSNP.MH.US.Outlier.thin2kb
less PBEfinalSNP.MH.US.Outlier.thin2kb.recode.vcf | grep -v '^#' | awk -v OFS="\t" '{start=$2-10000; end=$2+10000; print $1, start, end}' > MH.US.OutlierSNP.focalRegions.txt

vcftools --vcf landrace_palmarChico_PBEfinalSNP_v3.recode.vcf --positions MH.AN.commonOutlierSNP.txt --thin 2000 --recode --out PBEfinalSNP.MH.AN.Outlier.thin2kb
less PBEfinalSNP.MH.AN.Outlier.thin2kb.recode.vcf | grep -v '^#' | awk -v OFS="\t" '{start=$2-10000; end=$2+10000; print $1, start, end}' > MH.AN.OutlierSNP.focalRegions.txt

vcftools --vcf landrace_palmarChico_PBEfinalSNP_v3.recode.vcf --positions GH.AN.commonOutlierSNP.txt --thin 2000 --recode --out PBEfinalSNP.GH.AN.Outlier.thin2kb
less PBEfinalSNP.GH.AN.Outlier.thin2kb.recode.vcf | grep -v '^#' | awk -v OFS="\t" '{start=$2-10000; end=$2+10000; print $1, start, end}' > GH.AN.OutlierSNP.focalRegions.txt

vcftools --vcf landrace_palmarChico_PBEfinalSNP_v3.recode.vcf --positions US.AN.commonOutlierSNP.txt --thin 2000 --recode --out PBEfinalSNP.US.AN.Outlier.thin2kb
less PBEfinalSNP.US.AN.Outlier.thin2kb.recode.vcf | grep -v '^#' | awk -v OFS="\t" '{start=$2-10000; end=$2+10000; print $1, start, end}' > US.AN.OutlierSNP.focalRegions.txt

vcftools --vcf landrace_palmarChico_PBEfinalSNP_v3.recode.vcf --positions GH.US.commonOutlierSNP.txt --thin 2000 --recode --out PBEfinalSNP.GH.US.Outlier.thin2kb
less PBEfinalSNP.GH.US.Outlier.thin2kb.recode.vcf | grep -v '^#' | awk -v OFS="\t" '{start=$2-10000; end=$2+10000; print $1, start, end}' > GH.US.OutlierSNP.focalRegions.txt


cut -d ' ' -f1,2,5 MH.ML.PARV.GH.US.AN.SL.REFfrq  | awk -v FS=" " -v OFS="\t" '$2!=$3' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.GH.REFfrq
cut -f1,4 MH.ML.PARV.GH.US.AN.SL.REFfrq.neutral  | awk -v FS=" " -v OFS="\t" '$1!=$2' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.GH.REFfrq.neutral

cut -d ' ' -f1,2,7 MH.ML.PARV.GH.US.AN.SL.REFfrq  | awk -v FS=" " -v OFS="\t" '$2!=$3' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.AN.REFfrq
cut -f1,6 MH.ML.PARV.GH.US.AN.SL.REFfrq.neutral  | awk -v FS=" " -v OFS="\t" '$1!=$2' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.AN.REFfrq.neutral

cut -d ' ' -f1,2,6 MH.ML.PARV.GH.US.AN.SL.REFfrq  | awk -v FS=" " -v OFS="\t" '$2!=$3' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.US.REFfrq
cut -f1,5 MH.ML.PARV.GH.US.AN.SL.REFfrq.neutral  | awk -v FS=" " -v OFS="\t" '$1!=$2' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.US.REFfrq.neutral

cut -d ' ' -f1,5,6 MH.ML.PARV.GH.US.AN.SL.REFfrq  | awk -v FS=" " -v OFS="\t" '$2!=$3' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > GH.US.REFfrq
cut -f4,5 MH.ML.PARV.GH.US.AN.SL.REFfrq.neutral  | awk -v FS=" " -v OFS="\t" '$1!=$2' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > GH.US.REFfrq.neutral

cut -d ' ' -f1,5,7 MH.ML.PARV.GH.US.AN.SL.REFfrq  | awk -v FS=" " -v OFS="\t" '$2!=$3' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > GH.AN.REFfrq
cut -f4,6 MH.ML.PARV.GH.US.AN.SL.REFfrq.neutral  | awk -v FS=" " -v OFS="\t" '$1!=$2' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > GH.AN.REFfrq.neutral

cut -d ' ' -f1,6,7 MH.ML.PARV.GH.US.AN.SL.REFfrq  | awk -v FS=" " -v OFS="\t" '$2!=$3' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > US.AN.REFfrq
cut -f5,6 MH.ML.PARV.GH.US.AN.SL.REFfrq.neutral  | awk -v FS=" " -v OFS="\t" '$1!=$2' | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > US.AN.REFfrq.neutral


join -1 1 -2 1 <(awk -v OFS="\t" '{print $1"_"$2, $3, $4}' MH.GH.REFfrq  | sort -k1,1) <(awk -v OFS="\t" '{print $1"_"$2}' MH.GH.commonNeutralSNP.txt | sort -k1,1) | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.GH.REFfrq.neutral
join -1 1 -2 1 <(awk -v OFS="\t" '{print $1"_"$2, $3, $4}' MH.US.REFfrq  | sort -k1,1) <(awk -v OFS="\t" '{print $1"_"$2}' MH.US.commonNeutralSNP.txt | sort -k1,1) | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.US.REFfrq.neutral
join -1 1 -2 1 <(awk -v OFS="\t" '{print $1"_"$2, $3, $4}' MH.AN.REFfrq  | sort -k1,1) <(awk -v OFS="\t" '{print $1"_"$2}' MH.AN.commonNeutralSNP.txt | sort -k1,1) | sed 's/\_/\t/g' | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > MH.AN.REFfrq.neutral


