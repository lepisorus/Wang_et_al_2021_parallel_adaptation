df <- read.delim("MH.GH.US.AN.REFfrq.neutral", header=F)
head(df)
df.slim <- df[sample(nrow(df), 500000, replace=F), ]
head(df.slim)
length(df.slim$V1)
df.t <- t(df.slim[, 3:6])
sampleSizes = c(8, 6, 12, 10)
neutralF_filename = "neutralF"
allFreqs <- df.t
source("calcNeutralF.R")
