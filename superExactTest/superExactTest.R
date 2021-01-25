#sed '1d' MexHigh.outlierPBE.txt | awk -v OFS="_" '{print $1, $2}' | sort -k1,1 > MH.outlierSNP.pos
#sed '1d' GuaHigh.outlierPBE.txt | awk -v OFS="_" '{print $1, $2}' | sort -k1,1 > GH.outlierSNP.pos
#sed '1d' SW_US.outlierPBE.txt | awk -v OFS="_" '{print $1, $2}' | sort -k1,1 > US.outlierSNP.pos
#sed '1d' Andes.outlierPBE.txt | awk -v OFS="_" '{print $1, $2}' | sort -k1,1 > AN.outlierSNP.pos

library("SuperExactTest")

##SNP
MexHigh <- read.delim("MH.outlierSNP.pos", header=F)
GuaHigh <- read.delim("GH.outlierSNP.pos", header=F)
SW_US <- read.delim("US.outlierSNP.pos", header=F)
Andes <- read.delim("AN.outlierSNP.pos", header=F)

outlierSNP <- list(MH=as.character(MexHigh$V1), GH=as.character(GuaHigh$V1), US=as.character(SW_US$V1), AN=as.character(Andes$V1))

total= 1567351 #define the total number of genes/SNPs
res=supertest(outlierSNP, n=total)

pdf("outlierSNPintersectionEnrich.pdf")
plot(res, Layout="landscape", degree=2:4, sort.by="size", x.pos=c(0.1, 0.9), y.pos=c(0.025,0.9))
dev.off()

write.table(summary(res)$Table, file="outlierSNPintersectionEnrich.txt", row.names=FALSE, quote=F, sep="\t")

help("plot.msets")

##Gene
MexHigh <- read.delim("MexHigh.Outlier.fullGene.txt", header=F)
GuaHigh <- read.delim("GuaHigh.Outlier.fullGene.txt", header=F)
SW_US <- read.delim("SW_US.Outlier.fullGene.txt", header=F)
Andes <- read.delim("Andes.Outlier.fullGene.txt", header=F)


MexHigh <- read.delim("MexHighOutlier.1PercOutlier.fullGene.txt", header=F)
GuaHigh <- read.delim("GuaHighOutlier.1PercOutlier.fullGene.txt", header=F)
SW_US <- read.delim("SW_USOutlier.1PercOutlier.fullGene.txt", header=F)
Andes <- read.delim("AndesOutlier.1PercOutlier.fullGene.txt", header=F)

outlierGene <- list(MH=as.character(MexHigh$V1), GH=as.character(GuaHigh$V1), US=as.character(SW_US$V1), AN=as.character(Andes$V1))
total= 29160 #finalSNPs.fullGene.txt
res=supertest(outlierGene, n=total)

pdf("outlierGeneIntersectionEnrich.pdf")
plot(res, Layout="landscape", degree=2:4, sort.by="size", x.pos=c(0.1, 0.9), y.pos=c(0.025,0.9))
dev.off()

write.table(summary(res)$Table, file="outlierGeneintersectionEnrich.txt", row.names=FALSE, quote=F, sep="\t")

##flowering time gene
df <- read.delim("commonOutlierGenesInFloweringTime.txt")
MH <- df$geneID[df$MH==1]
GH <- df$geneID[df$GH==1]
US <- df$geneID[df$US==1]
AN <- df$geneID[df$AN==1]

outlierFlGene <- list(MH=as.character(MH), GH=as.character(GH), US=as.character(US), AN=as.character(AN))

total=920

res=supertest(outlierFlGene, n=total)

pdf("outlierFLgeneIntersectionEnrich.pdf")
plot(res, Layout="landscape", degree=2:4, sort.by="size", x.pos=c(0.08, 0.9), y.pos=c(0.025,0.9))
dev.off()

write.table(summary(res)$Table, file="outlierFLgeneIntersectionEnrich.txt", row.names=FALSE, quote=F, sep="\t")

##Glycerolipid genes
df <- read.delim("commonOutlierGenesInGlycerolipidPathway.txt")
MH <- df$geneID[df$MH==1]
GH <- df$geneID[df$GH==1]
US <- df$geneID[df$US==1]
AN <- df$geneID[df$AN==1]

outlierGLYGene <- list(MH=as.character(MH), GH=as.character(GH), US=as.character(US), AN=as.character(AN))

total=210

res=supertest(outlierGLYGene, n=total)

pdf("outlierGLYCEROLIPIDgeneIntersectionEnrich.pdf")
plot(res, Layout="landscape", degree=2:4, sort.by="size", x.pos=c(0.08, 0.9), y.pos=c(0.025,0.9))
dev.off()

write.table(summary(res)$Table, file="outlierGLYCEROLIPIDgeneIntersectionEnrich.txt", row.names=FALSE, quote=F, sep="\t")


###navarro all flowering time candidate genes
df <- read.delim("romero_navarro_flowering_outlier.txt", header=F)

MH <- read.delim("../outlierGenes/MexHigh.Outlier.fullGene.txt", header=F)
df$MH <- ifelse(df$V1 %in% MH$V1, 1, 0)
head(df)
table(df$MH)

US <- read.delim("../outlierGenes/SW_US.Outlier.fullGene.txt", header=F)
df$US <- ifelse(df$V1 %in% US$V1, 1, 0)
head(df)
table(df$US)

GH <- read.delim("../outlierGenes/GuaHigh.Outlier.fullGene.txt", header=F)
df$GH <- ifelse(df$V1 %in% GH$V1, 1, 0)
head(df)
table(df$GH)

AN <- read.delim("../outlierGenes/Andes.Outlier.fullGene.txt", header=F)
df$AN <- ifelse(df$V1 %in% AN$V1, 1, 0)
head(df)
table(df$AN)

write.table(df, file="romero_navarro_flowering_outlier_each_pop.txt", quote=F, sep="\t", row.names=F)
names(df)[1] <- "geneID"
MH <- df$geneID[df$MH==1]
GH <- df$geneID[df$GH==1]
US <- df$geneID[df$US==1]
AN <- df$geneID[df$AN==1]

outlierGene <- list(MH=as.character(MH), GH=as.character(GH), US=as.character(US), AN=as.character(AN))
total= length(df$geneID) #finalSNPs.fullGene.txt
res=supertest(outlierGene, n=total)

pdf("outlierNavarroFloweringGeneIntersectionEnrich.pdf")
plot(res, Layout="landscape", degree=2:4, sort.by="size", x.pos=c(0.1, 0.9), y.pos=c(0.025,0.9))
dev.off()

write.table(summary(res)$Table, file="outlierNavarroFloweringGeneintersectionEnrich.txt", row.names=FALSE, quote=F, sep="\t")

#Dong 2012 flowering time genes
df <- read.delim("Dong2012GeneID.txt", header=F)

MH <- read.delim("../outlierGenes/MexHigh.Outlier.fullGene.txt", header=F)
df$MH <- ifelse(df$V1 %in% MH$V1, 1, 0)
head(df)
table(df$MH)

US <- read.delim("../outlierGenes/SW_US.Outlier.fullGene.txt", header=F)
df$US <- ifelse(df$V1 %in% US$V1, 1, 0)
head(df)
table(df$US)

GH <- read.delim("../outlierGenes/GuaHigh.Outlier.fullGene.txt", header=F)
df$GH <- ifelse(df$V1 %in% GH$V1, 1, 0)
head(df)
table(df$GH)

AN <- read.delim("../outlierGenes/Andes.Outlier.fullGene.txt", header=F)
df$AN <- ifelse(df$V1 %in% AN$V1, 1, 0)
head(df)
table(df$AN)

write.table(df, file="Dong_flowering_outlier_each_pop.txt", quote=F, sep="\t", row.names=F)
names(df)[1] <- "geneID"
MH <- df$geneID[df$MH==1]
GH <- df$geneID[df$GH==1]
US <- df$geneID[df$US==1]
AN <- df$geneID[df$AN==1]

outlierGene <- list(MH=as.character(MH), GH=as.character(GH), US=as.character(US), AN=as.character(AN))
total= length(df$geneID) #finalSNPs.fullGene.txt
res=supertest(outlierGene, n=total)

pdf("outlierDongFloweringGeneIntersectionEnrich.pdf")
plot(res, Layout="landscape", degree=2:4, sort.by="size", x.pos=c(0.1, 0.9), y.pos=c(0.025,0.9))
dev.off()

write.table(summary(res)$Table, file="outlierDongFloweringGeneintersectionEnrich.txt", row.names=FALSE, quote=F, sep="\t")



