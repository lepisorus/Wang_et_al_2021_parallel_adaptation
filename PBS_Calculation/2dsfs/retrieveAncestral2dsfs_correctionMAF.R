library(data.table)
library(plyr)
library(tidyr)
library(reshape2)

##
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
df.slim <- df[, c(1:3, 6:7, 12, 15:16)]
head(df.slim)

df.slim2 <- subset(df.slim, df.slim$freqALT_MexLow.x>=0.3 & df.slim$freqALT_MexLow.x <= 0.7 & df.slim$parv.ALT >= 0.3 & df.slim$parv.ALT<=0.7)
table(df.slim2$SNPtype)

##
MH.GH.outlier <- read.delim("../MexHigh.SW_US.outlier.SNPtype")
head(MH.GH.outlier)
df <- merge(MH.GH.outlier, parv.freq, by="ID")
length(df$ID)
names(df)
df.slim <- df[, c(1:3, 6:7, 12, 15:16)]
head(df.slim)
df.slim2 <- subset(df.slim, df.slim$freqALT_MexLow.x>=0.3 & df.slim$freqALT_MexLow.x <= 0.7 & df.slim$parv.ALT >= 0.3 & df.slim$parv.ALT<=0.7)
table(df.slim2$SNPtype)

##
MH.GH.outlier <- read.delim("../MexHigh.Andes.outlier.SNPtype")
head(MH.GH.outlier)
df <- merge(MH.GH.outlier, parv.freq, by="ID")
length(df$ID)
names(df)
df.slim <- df[, c(1:3, 6:7, 13:14, 16, 19:20)]
head(df.slim)
df.slim2 <- subset(df.slim, df.slim$freqALT_MexLow>=0.3 & df.slim$freqALT_MexLow <= 0.7 & df.slim$parv.ALT >= 0.3 & df.slim$parv.ALT<=0.7 & df.slim$freqALT_SA_Low >=0.3 & df.slim$freqALT_SA_Low<=0.7)
table(df.slim2$SNPtype)

##
MH.GH.outlier <- read.delim("../Andes.GuaHigh.outlier.SNPtype")
head(MH.GH.outlier)
df <- merge(MH.GH.outlier, parv.freq, by="ID")
length(df$ID)
names(df)
df.slim <- df[, c(1:3, 6:7, 13:14, 16, 19:20)]
head(df.slim)
df.slim2 <- subset(df.slim, df.slim$freqALT_MexLow>=0.3 & df.slim$freqALT_MexLow <= 0.7 & df.slim$parv.ALT >= 0.3 & df.slim$parv.ALT<=0.7 & df.slim$freqALT_SA_Low >=0.3 & df.slim$freqALT_SA_Low<=0.7)
table(df.slim2$SNPtype)

##
MH.GH.outlier <- read.delim("../Andes.SW_US.outlier.SNPtype")
head(MH.GH.outlier)
df <- merge(MH.GH.outlier, parv.freq, by="ID")
length(df$ID)
names(df)
df.slim <- df[, c(1:3, 6:7, 13:14, 16, 19:20)]
head(df.slim)
df.slim2 <- subset(df.slim, df.slim$freqALT_MexLow>=0.3 & df.slim$freqALT_MexLow <= 0.7 & df.slim$parv.ALT >= 0.3 & df.slim$parv.ALT<=0.7 & df.slim$freqALT_SA_Low >=0.3 & df.slim$freqALT_SA_Low<=0.7)
table(df.slim2$SNPtype)

##
MH.GH.outlier <- read.delim("../SW_US.GuaHigh.outlier.SNPtype")
head(MH.GH.outlier)
df <- merge(MH.GH.outlier, parv.freq, by="ID")
length(df$ID)
names(df)
df.slim <- df[, c(1:3, 6:7, 12, 15:16)]
head(df.slim)
df.slim2 <- subset(df.slim, df.slim$freqALT_MexLow.x>=0.3 & df.slim$freqALT_MexLow.x <= 0.7 & df.slim$parv.ALT >= 0.3 & df.slim$parv.ALT<=0.7)
table(df.slim2$SNPtype)
