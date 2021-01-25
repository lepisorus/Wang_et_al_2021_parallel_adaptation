library(cowplot)
library(wesanderson)


df <- read.delim("PIF3freq.txt")
head(df)
df$pop <- factor(df$pop, levels=unique(df$pop))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("PIF3.1 (3:48809146)")
ggsave("PIF3.1SNPfreq.pdf")

df <- read.delim("ZCN8freq.txt")
head(df)
df$pop <- factor(df$pop, levels=unique(df$pop))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("ZCN8 (8:123030815)")
ggsave("ZCN8SNPfreq.pdf")

df <- read.delim("LUXfreq.txt")
head(df)
df$pop <- factor(df$pop, levels=unique(df$pop))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("LUX (8:156261065)")
ggsave("LUX_SNPfreq.pdf")

df <- read.delim("vgt1_131578990_missense.txt")
df$pop <- factor(df$pop, levels=unique(df$pop))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("vgt1 (8:131578990)")
ggsave("vgt1_131578990_SNPfreq.pdf")

df <- read.delim("vgt1_131579463_synonymous.txt")
df$pop <- factor(df$pop, levels=unique(df$pop))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("vgt1 (8:131579463)")
ggsave("vgt1_131579463_SNPfreq.pdf")

df <- read.delim("vgt1_131580179_3UTR.txt")
df$pop <- factor(df$pop, levels=unique(df$pop))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("vgt1 (8:131580179)")
ggsave("vgt1_131580179_SNPfreq.pdf")

df <- read.delim("PIF3freq2.txt")
cols <- c("US" = wes_palette("Darjeeling1")[3], "MH" = wes_palette("Darjeeling1")[1], "ML" = wes_palette("Royal1")[1], "GH" = wes_palette("Royal2")[3], "SL"="black", "AN"="purple", "PARV"="lightGreen")
head(df)
df$pop <- factor(df$pop, levels=c("US", "MH", "ML", "GH", "SL", "AN", "PARV"))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("PIF3.1 (3:48809139)")+
scale_fill_manual(values=cols)+
ylab("frequency")
ggsave("PIF3.1SNPfreq2.pdf")

df <- read.delim("PIF3freq.txt")
cols <- c("US" = wes_palette("Darjeeling1")[3], "MH" = wes_palette("Darjeeling1")[1], "ML" = wes_palette("Royal1")[1], "GH" = wes_palette("Royal2")[3], "SL"="black", "AN"="purple", "PARV"="lightGreen")
head(df)
df$pop <- factor(df$pop, levels=c("US", "MH", "ML", "GH", "SL", "AN", "PARV"))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + ggtitle("PIF3.1 (3:48809139)")+
scale_fill_manual(values=cols)+
ylab("frequency")+
theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), legend.position="none", plot.title = element_text(size=20))
ggsave("PIF3.1SNPfreq.pdf")

df <- read.delim("PSL1freq.txt")
cols <- c("US" = wes_palette("Darjeeling1")[3], "MH" = wes_palette("Darjeeling1")[1], "ML" = wes_palette("Royal1")[1], "GH" = wes_palette("Royal2")[3], "SL"="black", "AN"="purple", "PARV"="lightGreen")
head(df)
df$pop <- factor(df$pop, levels=c("US", "MH", "ML", "GH", "SL", "AN", "PARV"))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + ggtitle("PLS1 (2:230277873)")+
scale_fill_manual(values=cols)+
ylab("frequency") +
theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), legend.position="none", plot.title = element_text(size=20))
ggsave("PLS1SNPfreq.pdf")