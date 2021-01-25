library(data.table)
library(plyr)
library(tidyr)
library(reshape2)

##read in the common neutral SNPs
neutral <- read.delim("MexHigh.SW_US.neutralPBE.txt")
head(neutral)
neutral["ID"] <- paste(neutral$CHROM, neutral$POS, sep="_")

## find out the neutral SNP type (co, anti, no-dire)
SNPtype <- read.delim("../MexHigh.SW_US.parallel.SNPtype")
head(SNPtype)
names(SNPtype)
SNPtype.slim <- SNPtype[, c(1,7,12)]
head(SNPtype.slim)
neutral.SNPtype <- merge(neutral, SNPtype.slim, by="ID")
length(neutral$ID)
length(neutral.SNPtype$ID)
head(neutral.SNPtype)

## combine parv frequency into the big file
parv.freq <- fread("parv.slim.frq")
names(parv.freq)[3:4] <- c("parv.REF", "parv.ALT")
parv.freq <- as.data.frame(parv.freq)
parv.freq["ID"] <- paste(parv.freq$V1, parv.freq$V2, sep="_")
head(parv.freq)
neutral.SNPtype2 <- merge(neutral.SNPtype, parv.freq, by="ID")
length(neutral.SNPtype2$ID)
names(neutral.SNPtype2)
neutral.SNPtype3 <- neutral.SNPtype2[, c(1,4,9,5)]
head(neutral.SNPtype3)

# bin the ancestral 2d sfs in neutral SNPs
neutral.SNPtype3$MexLow.ALT <- cut(neutral.SNPtype3$freqALT_MexLow.x, breaks=seq(0, 1, 0.1), labels=seq(0, 0.9, 0.1))
neutral.SNPtype3$parv.ALT.new <- cut(neutral.SNPtype3$parv.ALT, breaks=seq(0, 1, 0.1), labels=seq(0, 0.9, 0.1))

neutral.SNPtype3$parv.ALT.new <- as.numeric(as.character(neutral.SNPtype3$parv.ALT.new))
neutral.SNPtype3$parv.ALT.new <- ifelse(neutral.SNPtype3$parv.ALT==0, 0, neutral.SNPtype3$parv.ALT.new)

neutral.SNPtype3$MexLow.ALT <- as.numeric(as.character(neutral.SNPtype3$MexLow.ALT))
neutral.SNPtype3$MexLow.ALT <- ifelse(neutral.SNPtype3$freqALT_MexLow.x==0, 0, neutral.SNPtype3$MexLow.ALT)

head(neutral.SNPtype3)
write.table(neutral.SNPtype3[, c(1,5,6,4)], file="MH_SW_US.neutralSNP.ancestral2dsfs.txt", quote=F, sep="\t", row.names=F)


# find out the 2dsfs in outlier SNPs
MH.GH.outlier <- read.delim("../MexHigh.SW_US.outlier.SNPtype")
head(MH.GH.outlier)
parv.freq <- fread("parv.slim.frq")
head(parv.freq)
names(parv.freq)[3:4] <- c("parv.REF", "parv.ALT")
parv.freq <- as.data.frame(parv.freq)
parv.freq["ID"] <- paste(parv.freq$V1, parv.freq$V2, sep="_")
head(parv.freq)
head(MH.GH.outlier)
length(MH.GH.outlier$ID)
df <- merge(MH.GH.outlier, parv.freq, by="ID")
length(df$ID)
names(df)
df.slim <- df[, c(1:3, 6:7, 15:16)]
head(df.slim)
df.slim$MexLow.ALT <- cut(df.slim$freqALT_MexLow.x, breaks=seq(0, 1, 0.1), labels=seq(0, 0.9, 0.1))
df.slim$parv.ALT.new <- cut(df.slim$parv.ALT, breaks=seq(0, 1, 0.1), labels=seq(0, 0.9, 0.1))
df.slim$parv.ALT.new <- as.numeric(as.character(df.slim$parv.ALT.new))
df.slim$parv.ALT.new <- ifelse(df.slim$parv.ALT==0, 0, df.slim$parv.ALT.new)
df.slim$MexLow.ALT <- as.numeric(as.character(df.slim$MexLow.ALT))
df.slim$MexLow.ALT <- ifelse(df.slim$freqALT_MexLow.x==0, 0, df.slim$MexLow.ALT)
head(df.slim)
summary(df.slim$MexLow.ALT)
summary(df.slim$parv.ALT.new)
names(df.slim)
df.slim2 <- df.slim[, c(1, 8:9)]
t <- count(df.slim2, c('MexLow.ALT', 'parv.ALT.new'))
write.table(t, file="MH.SW_US.outlier.ancestral2dsfs.txt", quote=F, sep="\t", row.names=F)

## find out neutral SNPs with the same ancestral 2dsfs as in the outlier set
t["2dsfs"] <- paste(t$MexLow.ALT, t$parv.ALT.new, sep="_")
head(t)
neutral.SNPtype3["2dsfs"] <- paste(neutral.SNPtype3$MexLow.ALT, neutral.SNPtype3$parv.ALT.new, sep="_")
head(neutral.SNPtype3)
test <- merge(t, neutral.SNPtype3, by="2dsfs")

head(test)

#test.slim <- test2[,c(1,4:5,8)]
test.slim <- test[,c(1,4:5,8)]
head(test.slim)
names(test.slim)[1] <- "sfs"
length(test.slim$freq)

## split the data frame into list according to 2dsfs
test.list <- split(test.slim, test.slim$sfs)
length(test.list)

# for each allele frequency type, sample the same amount of neutral SNPs as that in outlier set
test.list2 <- sapply(test.list, function(x)sample(x$ID, unique(x$freq), replace=TRUE))
df <- ldply (test.list2, data.frame)
head(df)
names(df) <- c("sfs", "ID")

# got the SNPtype of the randomly sampled IDs ---> the same 2dsfs in neutral set and outlier set
head(df)
df.new <- merge(df, test.slim, by="ID")
head()
summary(df.new$SNPtype)

