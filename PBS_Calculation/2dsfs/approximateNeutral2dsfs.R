library(data.table)
library(plyr)
library(tidyr)
library(reshape2)

neutral <- read.delim("MexHigh.GuaHigh.neutralPBE.txt")
head(neutral)
neutral["ID"] <- paste(neutral$CHROM, neutral$POS, sep="_")

SNPtype <- read.delim("MexHigh.GuaHigh.parallel.SNPtype")
head(SNPtype)
names(SNPtype)
SNPtype.slim <- SNPtype[, c(1,7,12)]
head(SNPtype.slim)
neutral.SNPtype <- merge(neutral, SNPtype.slim, by="ID")
length(neutral$ID)
length(neutral.SNPtype$ID)
head(neutral.SNPtype)

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


neutral.SNPtype3$MexLow.ALT <- cut(neutral.SNPtype3$freqALT_MexLow.x, breaks=seq(0, 1, 0.1), labels=seq(0, 0.9, 0.1))
neutral.SNPtype3$parv.ALT.new <- cut(neutral.SNPtype3$parv.ALT, breaks=seq(0, 1, 0.1), labels=seq(0, 0.9, 0.1))

neutral.SNPtype3$parv.ALT.new <- as.numeric(as.character(neutral.SNPtype3$parv.ALT.new))
neutral.SNPtype3$parv.ALT.new <- ifelse(neutral.SNPtype3$parv.ALT==0, 0, neutral.SNPtype3$parv.ALT.new)

neutral.SNPtype3$MexLow.ALT <- as.numeric(as.character(neutral.SNPtype3$MexLow.ALT))
neutral.SNPtype3$MexLow.ALT <- ifelse(neutral.SNPtype3$freqALT_MexLow.x==0, 0, neutral.SNPtype3$MexLow.ALT)

head(neutral.SNPtype3)
write.table(neutral.SNPtype3[, c(1,5,6,4)], file="neutralSNP.ancestral2dsfs.txt", quote=F, sep="\t", row.names=F)


