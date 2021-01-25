library(data.table)
library(optparse)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
make_option(c('-o','--out_file1'), action='store', type='character', default=NULL, help='Output file'),
make_option(c('-s','--out_file2'), action='store', type='character', default=NULL, help='Output file'),
make_option(c('-n','--out_file3'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

df <- fread(opt$in_file)
head(df)
names(df)[3:5] <- c("Fst01", "Fst02", "Fst12")
head(df)
df <- as.data.frame(df)

## assign any Fst values < 0 as 0
df$Fst01 <- ifelse(df$Fst01 < 0, 0, df$Fst01)
df$Fst02 <- ifelse(df$Fst02 < 0, 0, df$Fst02)
df$Fst12 <- ifelse(df$Fst12 < 0, 0, df$Fst12)

df["T01"] <- -log(1-df$Fst01)
df["T02"] <- -log(1-df$Fst02)
df["T12"] <- -log(1-df$Fst12)
head(df)
df["PBS0"] <- (df$T01+df$T02-df$T12)/2
head(df)
df.new <- df[is.finite(df$T12) & is.finite(df$PBS0), ]
summary(df.new$PBS0)
summary(df.new$T12)

df.new1 <- df[is.finite(df$T12),  ]
T12.median <- median(df.new1$T12)

df.new2 <- df[is.finite(df$PBS0),  ]
PBS.median <- median(df.new$PBS0)

df.new["PBE0"] <- df.new$PBS0 - (df.new$T12*(PBS.median/T12.median))
head(df.new)
summary(df.new$PBE0)


write.table(df.new, file=opt$out_file1, quote=F, sep="\t", row.names=F)


df.outlier <- subset(df.new[, c(1,2)], df.new$PBE0 > quantile(df.new$PBE0, 0.95))
df.neutral <- subset(df.new[, c(1,2)], df.new$PBE0 < quantile(df.new$PBE0, 0.75) & df.new$PBE0 > quantile(df.new$PBE0, 0.25))

write.table(df.outlier, file=opt$out_file2, quote=F, sep="\t", row.names=F)
write.table(df.neutral, file=opt$out_file3, quote=F, sep="\t", row.names=F)

