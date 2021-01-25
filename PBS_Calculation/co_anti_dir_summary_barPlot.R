library(cowplot)
library(reshape2)

df <- read.delim("co_anti_dir_summary.txt")
head(df)
df.slim <- df[, c(1,2,7,8)]
head(df.slim)
df.long <- melt(df.slim, id.var=c("pops", "type"))
df.long$SNPtype <- ifelse(df.long$variable=="coDir.freq", "co-direction", "anti-direction")
head(df.long)
names(df.long)[4] <- "frequency"
table(df.long$pops)
head(df.long)

df.long["types"] <- ifelse(df.long$type=="outlier", "o", "n")
df.long["populations"] <- ifelse(df.long$pops=="GuaHigh_Andes", "GH_AD", "unknown")
df.long$populations <- ifelse(df.long$pops=="GuaHigh_SW_US", "GH_SU", df.long$populations)
df.long$populations <- ifelse(df.long$pops=="MexHigh_Andes", "MH_AD", df.long$populations)
df.long$populations <- ifelse(df.long$pops=="MexHigh_GuaHigh", "MH_GH", df.long$populations)
df.long$populations <- ifelse(df.long$pops=="MexHigh_SW_US", "MH_SU", df.long$populations)
df.long$populations <- ifelse(df.long$pops=="SW_US_Andes", "SU_AD", df.long$populations)

p <- ggplot(df.long, aes(x=types, y=frequency, fill=SNPtype)) + geom_bar(stat="identity", width=0.5) + 
facet_wrap(~populations,nrow = 1)  + 
#scale_color_manual(values=c("skyblue", "orange"))+
xlab("") + ylab("frequency")
p
ggsave("co_anti_dir_summary.pdf")

