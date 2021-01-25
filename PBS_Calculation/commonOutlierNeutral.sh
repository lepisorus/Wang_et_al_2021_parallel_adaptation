cat MexHigh.outlierPBE.txt GuaHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > MexHigh.GuaHigh.common.outlierPBE.txt
cat MexHigh.neutralPBE.txt GuaHigh.neutralPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > MexHigh.GuaHigh.neutralPBE.txt

cat MexHigh.outlierPBE.txt SW_US.outlierPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > MexHigh.SW_US.common.outlierPBE.txt
cat MexHigh.neutralPBE.txt SW_US.neutralPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > MexHigh.SW_US.neutralPBE.txt

cat MexHigh.outlierPBE.txt Andes.outlierPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > MexHigh.Andes.common.outlierPBE.txt
cat MexHigh.neutralPBE.txt Andes.neutralPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > MexHigh.Andes.neutralPBE.txt

cat SW_US.outlierPBE.txt GuaHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > SW_US.GuaHigh.common.outlierPBE.txt
cat SW_US.neutralPBE.txt GuaHigh.neutralPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > SW_US.GuaHigh.neutralPBE.txt

cat Andes.outlierPBE.txt GuaHigh.outlierPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > Andes.GuaHigh.common.outlierPBE.txt
cat Andes.neutralPBE.txt GuaHigh.neutralPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > Andes.GuaHigh.neutralPBE.txt

cat Andes.outlierPBE.txt SW_US.outlierPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > Andes.SW_US.common.outlierPBE.txt
cat Andes.neutralPBE.txt SW_US.neutralPBE.txt | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==2{print $2, $3}' > Andes.SW_US.neutralPBE.txt
