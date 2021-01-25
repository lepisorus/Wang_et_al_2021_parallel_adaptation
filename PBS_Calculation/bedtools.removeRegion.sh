bedtools intersect -v -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2, $3, $4, $5}' MexHigh.weir.fst) -b MexHigh.filterRegion.txt | cut -f1,2,4-6 > MexHigh.weir.fst.rmReg
bedtools intersect -v -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2, $3, $4, $5}' GuaHigh.weir.fst) -b GuaHigh.filterRegion.txt | cut -f1,2,4-6 > GuaHigh.weir.fst.rmReg
bedtools intersect -v -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2, $3, $4, $5}' SW_US.weir.fst) -b SW_US.filterRegion.txt | cut -f1,2,4-6 > SW_US.weir.fst.rmReg
bedtools intersect -v -a <(awk -v OFS="\t" 'NR!=1{print $1, $2, $2, $3, $4, $5}' Andes.weir.fst) -b Andes.filterRegion.txt | cut -f1,2,4-6 > Andes.weir.fst.rmReg
