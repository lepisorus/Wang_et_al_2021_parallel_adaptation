gene <- read.delim("Dong2012GeneID.txt", header=F)
head(gene)
within <- read.delim("finalSNPs.withinGene.txt", header=F)
head(within)
within1 <- subset(within, within$V7 %in% gene$V1)

outside <- read.delim("finalSNPs.10kbGene.txt", header=F)
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

df <- read.delim("../MexHigh.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("FLoutlierPBEdis_MexHigh.pdf") ##p=0.001
hist(PBE.dist, col="black", breaks=50, xlim=c(0.016, 0.064), xlab="PBE", main="")
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off()

#GuaHigh

df <- read.delim("../GuaHigh.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("FLoutlierPBEdis_GuaHigh.pdf") ##0.0006
hist(PBE.dist, col="black", breaks=50, xlim=c(0.029, 0.076), xlab="PBE", main="")
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off()

#SW_US

df <- read.delim("../SW_US.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("FLoutlierPBEdis_SW_US.pdf") ##0.021
hist(PBE.dist, col="black", breaks=50, xlim=c(0.014, 0.059), xlab="PBE", main="")
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off()

##

df <- read.delim("../Andes.allPBE.txt")
head(df)
df["ID"] <- paste(df$CHROM, df$POS, sep="_")
df.FLoutlier <- merge(df, FLoutlier, by="ID")
length(df.FLoutlier$ID)
PBE.dist <- replicate(10000, mean(sample(df$PBE0, length(df.FLoutlier$PBE0), replace=F)))
summary(df.FLoutlier$PBE0)
summary(PBE.dist)

pdf("FLoutlierPBEdis_Andes.pdf")
hist(PBE.dist, col="black", breaks=50, xlim=c(0.034, 0.17), xlab="PBE", main="")
abline(v=mean(df.FLoutlier$PBE), col="red", lwd=2)
dev.off() #p < 0.0001

