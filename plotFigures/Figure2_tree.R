library(data.table)
library(cowplot)

df <- fread("../MexHigh.allPBE.txt")
df <- as.data.frame(df)
df["SNPname"] <- paste(df$CHROM, df$POS, sep="_")
cutoff <- quantile(df$PBE0, 0.95)
cutoff

PIF3.1 <- subset(df, df$CHROM==3 & df$POS>=48000000 & df$POS<=50000000)
#neutral point 49964925 
pdf("MH_PIF3.1_PBEdis.pdf")
ggplot(data=PIF3.1, aes(x=POS/1000000, y=PBE0)) +  xlab("physical position (mb)") + ylab("PBE") + ggtitle("chr3 PIF3.1") + geom_hline(yintercept=0.4, color="red")+
geom_rect(aes(xmin=48.806991, xmax=48.810673, ymin=-Inf, ymax=Inf), color="lightgrey", alpha=0.1)+
geom_point(size=1) 
dev.off()


library(phangorn)

df <- matrix(c(1.95, 0, 0), nrow=3)
d <- dist(df, diag=F)
tree <- nj(d)
pdf("PIF3.1tree.pdf")
plot(tree)
dev.off()

nnet <- neighborNet(d)
pdf("PIF3.1tree2.pdf")
plot(nnet, "2D")
dev.off()

#
df2 <- matrix(c(-0.387395184800148, 0.569717141594185, 0.569717141594185), nrow=3)
d <- dist(df2, diag=F)
tree <- nj(d)
pdf("negativePBEtree.pdf")
plot(tree)
dev.off()

#
df3 <- matrix(c(0.0911609783970173, -0.0911609783970175, 0.0911609783970175), nrow=3)
d <- dist(df3, diag=F)
tree <- nj(d)
pdf("neutralPBEtree.pdf")
plot(tree)
dev.off()

library(cowplot)
df <- read.delim("PIF3freq2.txt")
head(df)
df$pop <- factor(df$pop, levels=unique(df$pop))
ggplot(df, aes(x=pop, y=freq, fill=pop)) + geom_bar(stat="identity", width=0.5) + 
scale_x_discrete(breaks=df$pop[nchar(as.character(df$pop))!=1]) +
xlab("") + theme(legend.position="none") + ggtitle("PIF3.1 (3:48809139)") +
ylab("frequency")
ggsave("PIF3.1SNPfreq_2.pdf")