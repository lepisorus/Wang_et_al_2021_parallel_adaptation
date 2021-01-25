library(ggplot2)

#MH VS AN
df <- read.delim("ML_AN_MH.fd.modelSummary.txt", header=F)
head(df)
header <- read.delim("header.txt", header=F)
header
header <- as.vector(t(header))
header
names(df) <- header
head(df)
df2 <- df[complete.cases(df), ]
head(df2)
length(df2$model)
length(df$model)

mig <- subset(df2, df2$model=="mig")
other <- subset(df2, df2$model!="mig")
wilcox.test(mig$fdM, other$fdM, alternative="greater")
df2 <-  subset(df2, df2$model!="neutral")
ggplot(data=df2, aes(x=model, y=fdM, fill=model, color=model)) + geom_violin() + ggtitle("MH vs AN convergent SNPs")+
stat_summary(fun.y=mean, geom="point", size=2, color="black")+
theme(legend.position="none")
ggsave("MH.AN.fdM.model.violin.pdf")

#GH VS AN
df <- read.delim("ML_AN_GH.fd.modelSummary.txt", header=F)
head(df)
header <- read.delim("header.txt", header=F)
header
header <- as.vector(t(header))
header
names(df) <- header
head(df)
df2 <- df[complete.cases(df), ]
head(df2)
length(df2$model)
length(df$model)

mig <- subset(df2, df2$model=="mig")
other <- subset(df2, df2$model!="mig")
wilcox.test(mig$fdM, other$fdM, alternative="greater")
df2 <-  subset(df2, df2$model!="neutral")
ggplot(data=df2, aes(x=model, y=fdM, fill=model, color=model)) + geom_violin() + ggtitle("GH vs AN convergent SNPs")+
stat_summary(fun.y=mean, geom="point", size=2, color="black")+
theme(legend.position="none")
ggsave("GH.AN.fdM.model.violin.pdf")

#US VS AN
df <- read.delim("ML_AN_US.fd.modelSummary.txt", header=F)
head(df)
header <- read.delim("header.txt", header=F)
header
header <- as.vector(t(header))
header
names(df) <- header
head(df)
df2 <- df[complete.cases(df), ]
head(df2)
length(df2$model)
length(df$model)

mig <- subset(df2, df2$model=="mig")
other <- subset(df2, df2$model!="mig")
wilcox.test(mig$fdM, other$fdM, alternative="greater")
df2 <-  subset(df2, df2$model!="neutral")
ggplot(data=df2, aes(x=model, y=fdM, fill=model, color=model)) + geom_violin() + ggtitle("US vs AN convergent SNPs")+
stat_summary(fun.y=mean, geom="point", size=2, color="black")+
theme(legend.position="none")
ggsave("US.AN.fdM.model.violin.pdf")
