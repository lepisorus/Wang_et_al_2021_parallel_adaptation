gene <- read.delim("glycerolipidGenes.txt", header=F)
head(gene)
within <- read.delim("../finalSNPs.withinGene.txt", header=F)
head(within)
within1 <- subset(within, within$V7 %in% gene$V1)

outside <- read.delim("../finalSNPs.10kbGene.txt", header=F)
head(outside)
outside1 <- subset(outside, outside$V7 %in% gene$V1)
length(outside1$V1)
length(within1$V1)
ID1 <- paste(within1$V1, within1$V2, sep="_")
head(ID1)
ID2 <- paste(outside1$V1, outside1$V2, sep="_")
head(ID2)
ID <- c(ID1, ID2)
FLoutlier <- as.data.frame(ID)

df <- read.delim("../../MexHigh.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("GLYoutlierPBEdis_MexHigh.pdf") ##p=0.001
hist(PBE.dist, col="black", breaks=50, xlim=c(0.025, 0.072), xlab="PBE", main="MH p<0.0001", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off()

#GuaHigh

df <- read.delim("../../GuaHigh.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("GLYoutlierPBEdis_GuaHigh.pdf") ##0.0006
hist(PBE.dist, col="black", breaks=50, xlim=c(0.04, 0.07), xlab="PBE", main="GH p<0.0001", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off()

#SW_US

df <- read.delim("../../SW_US.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("GLYoutlierPBEdis_SW_US.pdf") ##0.021
hist(PBE.dist, col="black", breaks=50, xlim=c(0.024, 0.062), xlab="PBE", main="US p<0.0001", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off()

##

df <- read.delim("../../Andes.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("GLYoutlierPBEdis_Andes.pdf")
hist(PBE.dist, col="black", breaks=50, xlim=c(0.056, 0.096), xlab="PBE", main="AN p<0.0001", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off() #p < 0.0001

