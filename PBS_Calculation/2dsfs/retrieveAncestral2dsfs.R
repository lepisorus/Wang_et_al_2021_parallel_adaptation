library(data.table)
library(plyr)
library(tidyr)
library(reshape2)

MH.GH.outlier <- read.delim("../MexHigh.GuaHigh.outlier.SNPtype")
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


names(df.slim)
df.slim2 <- df.slim[, c(1, 8:9)]

t <- count(df.slim2, c('MexLow.ALT', 'parv.ALT.new'))

write.table(t, file="MH.GH.outlier.ancestral2dsfs.txt", quote=F, sep="\t", row.names=F)




