library("SuperExactTest")
##Gene
MexHigh <- read.delim("MexHigh.outlierPBE.prune.txt2", header=T)
GuaHigh <- read.delim("GuaHigh.outlierPBE.prune.txt2", header=T)
SW_US <- read.delim("SW_US.outlierPBE.prune.txt2", header=T)
Andes <- read.delim("Andes.outlierPBE.prune.txt2", header=T)

outlierGene <- list(MH=as.character(MexHigh$CHROM_POS), GH=as.character(GuaHigh$CHROM_POS), US=as.character(SW_US$CHROM_POS), AN=as.character(Andes$CHROM_POS))
total= 515422 #finalSNPs.fullGene.txt
res=supertest(outlierGene, n=total)

pdf("outlierSNPIntersectionEnrich.pdf")
plot(res, Layout="landscape", degree=2:4, sort.by="size", x.pos=c(0.1, 0.9), y.pos=c(0.025,0.9), margin=c(1.5,5,1.5,5))
dev.off()

write.table(summary(res)$Table, file="outlierSNPintersectionEnrich.txt", row.names=FALSE, quote=F, sep="\t")
