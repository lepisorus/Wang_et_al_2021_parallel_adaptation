library(dgconstraint)

##for all selected genes
#this file contains genes involved in the final SNP dataset
df <- read.delim("finalSNPs.fullGene.txt", header=F)
head(df)
MexHigh <- read.delim("MexHigh.Outlier.fullGene.txt", header=F)
GuaHigh <- read.delim("GuaHigh.Outlier.fullGene.txt", header=F)
SW_US <- read.delim("SW_US.Outlier.fullGene.txt", header=F)
Andes <- read.delim("Andes.Outlier.fullGene.txt", header=F)
df["MH"] <- ifelse(df$V1 %in% MexHigh$V1, "1", "0")
df["GH"] <- ifelse(df$V1 %in% GuaHigh$V1, "1", "0")
df["US"] <- ifelse(df$V1 %in% SW_US$V1, "1", "0")
df["AN"] <- ifelse(df$V1 %in% Andes$V1, "1", "0")

##doesnot like the mean value of C_hyper returned, changed the function a bit to return the matrix of C_hyper
pairwise_c_hyper <- function (input, na.rm = F){
  numcol <- ncol (input)
  results_c_hyper <- array (NA, c (numcol,numcol))
  for(loop1 in 1:(numcol - 1)){
    for(loop2 in loop1:numcol){
      if (loop1 != loop2){
        ax <- sum (input[,loop1], na.rm = na.rm)
        ay <- sum (input[,loop2], na.rm = na.rm)
        g0 <- nrow (input)
        sd_hyp <- sqrt((ax*ay)*(g0-ax)*(g0-ay)/(g0^2*(g0-1)))
        exp_hyp <- ax * ay / g0
        obs_hyp <- sum (input[,loop1] == 1 & input[,loop2] == 1, na.rm = na.rm)
        if (sd_hyp != 0){
          results_c_hyper[loop1,loop2] <- (obs_hyp - exp_hyp) / sd_hyp
        } else {
          results_c_hyper[loop1,loop2] <- 0
          warning ('Some pairwise contrasts have no shared adapted loci')
        }
      }
    }
  }
results_c_hyper
}


pairwise_c_chisq <- function (input,num_permute = 10000,na.rm = F){

  numcol <- ncol (input)
  c_score <- array(NA, c(numcol,numcol))

  for (loop1 in 1:(numcol-1)){
    for (loop2 in (loop1+1):numcol){
      results_chisq <- array (NA,num_permute)
      input_sub <- as.matrix (input[,c(loop1,loop2)])
      obs1 <- rowSums (input_sub,na.rm = na.rm)
      exp1 <- array (mean(obs1), length(obs1))
      chisq1 <- (obs1 - exp1)^2 / exp1

      for (loop3 in 1:num_permute){
        input_sub[,1] <- sample (input_sub[,1],nrow (input_sub),replace = F)
        input_sub[,2] <- sample (input_sub[,2],nrow (input_sub),replace = F)

        obs2 <- rowSums(input_sub,na.rm = na.rm)
        exp2 <- array (mean(obs2),length(obs2))
        chisq2 <- (obs2 - exp2)^2/exp2

        results_chisq[loop3] <- sum (chisq2)
      }
      c_score [loop1,loop2]  <- (sum (chisq1) - mean (results_chisq)) / sd (results_chisq)
    }
  }

 c_score

}

row.names(df) <- df$V1
head(df)
df.new <- data.matrix(df[, 2:5])
write.table(df.new, file="allGene.adaptiveMatrix.txt", quote=F, sep="\t", row.names=F)

df.new <- read.table("allGene.adaptiveMatrix.txt", header=T)

MH_GH <- df.new[, 1:2]
MH_US <- df.new[, c(1,3)]
MH_AN <- df.new[, c(1,4)]
GH_US <- df.new[, 2:3]
GH_AN <- df.new[, c(2,4)]
US_AN <- df.new[, 3:4]

#pairwise C_chisq
C_chisq_MH_GH <- pairwise_c_chisq(MH_GH)
C_chisq_MH_US <- pairwise_c_chisq(MH_US)
C_chisq_MH_AN <- pairwise_c_chisq(MH_AN)
C_chisq_GH_US <- pairwise_c_chisq(GH_US)
C_chisq_GH_AN <- pairwise_c_chisq(GH_AN)
C_chisq_US_AN <- pairwise_c_chisq(US_AN)

C_chisq <- as.data.frame(cbind(C_chisq_MH_GH, C_chisq_MH_US, C_chisq_MH_AN, C_chisq_GH_US, C_chisq_GH_AN, C_chisq_US_AN))

write.table(C_chisq, file="C_chisq.txt", quote=F, sep="\t", row.names=F)


#pairwise C_chisq_p
C_chisq_p_MH_GH <- allwise_p_chisq(MH_GH)
C_chisq_p_MH_US <- allwise_p_chisq(MH_US)
C_chisq_p_MH_AN <- allwise_p_chisq(MH_AN)
C_chisq_p_GH_US <- allwise_p_chisq(GH_US)
C_chisq_p_GH_AN <- allwise_p_chisq(GH_AN)
C_chisq_p_US_AN <- allwise_p_chisq(US_AN)

C_chisq_p <- as.data.frame(cbind(C_chisq_p_MH_GH, C_chisq_p_MH_US, C_chisq_p_MH_AN, C_chisq_p_GH_US, C_chisq_p_GH_AN, C_chisq_p_US_AN))

write.table(C_chisq_p, file="C_chisq_p.txt", quote=F, sep="\t", row.names=F)

##null hypothesis: low GP redundancy:

df.new["sum"] <- as.numeric(df.new$MH) + as.numeric(df.new$GH) + as.numeric(df.new$US) +as.numeric(df.new$AN)
df2 <- subset(df.new, df.new$sum != 0)
head(df2)
df.slim <- data.matrix(df2)

MH_GH <- df.slim[, 1:2]
MH_US <- df.slim[, c(1,3)]
MH_AN <- df.slim[, c(1,4)]
GH_US <- df.slim[, 2:3]
GH_AN <- df.slim[, c(2,4)]
US_AN <- df.slim[, 3:4]

C_chisq_MH_GH <- pairwise_c_chisq(MH_GH)
C_chisq_MH_US <- pairwise_c_chisq(MH_US)
C_chisq_MH_AN <- pairwise_c_chisq(MH_AN)
C_chisq_GH_US <- pairwise_c_chisq(GH_US)
C_chisq_GH_AN <- pairwise_c_chisq(GH_AN)
C_chisq_US_AN <- pairwise_c_chisq(US_AN)

C_chisq <- as.data.frame(cbind(C_chisq_MH_GH, C_chisq_MH_US, C_chisq_MH_AN, C_chisq_GH_US, C_chisq_GH_AN, C_chisq_US_AN))

write.table(C_chisq, file="C_chisq_slim.txt", quote=F, sep="\t", row.names=F)


C_chisq_slim_p_MH_GH <- allwise_p_chisq(MH_GH)
C_chisq_slim_p_MH_US <- allwise_p_chisq(MH_US)
C_chisq_slim_p_MH_AN <- allwise_p_chisq(MH_AN)
C_chisq_slim_p_GH_US <- allwise_p_chisq(GH_US)
C_chisq_slim_p_GH_AN <- allwise_p_chisq(GH_AN)
C_chisq_slim_p_US_AN <- allwise_p_chisq(US_AN)

C_chisq_slim_p <- as.data.frame(cbind(C_chisq_slim_p_MH_GH, C_chisq_slim_p_MH_US, C_chisq_slim_p_MH_AN, C_chisq_slim_p_GH_US, C_chisq_slim_p_GH_AN, C_chisq_slim_p_US_AN))

write.table(C_chisq_slim_p, file="C_chisq_slim_p.txt", quote=F, sep="\t", row.names=F)


####### flowering time  #######
df <- read.delim("commonOutlierGenesInFloweringTime.txt", header=T)
row.names(df) <- df$V1
head(df)
df.new <- data.matrix(df[, 2:5])


C_chisq <- pairwise_c_chisq(df.new)
row.names(C_chisq) <- c("MH", "GH", "US", "AN")
write.table(C_chisq, file="pairwise_C_chisq_flowering.txt", quote=F, sep="\t", row.names=F)

MH_GH <- df.new[, 1:2]
MH_US <- df.new[, c(1,3)]
MH_AN <- df.new[, c(1,4)]
GH_US <- df.new[, 2:3]
GH_AN <- df.new[, c(2,4)]
US_AN <- df.new[, 3:4]

#pairwise C_chisq
C_chisq_MH_GH <- pairwise_c_chisq(MH_GH)
C_chisq_MH_US <- pairwise_c_chisq(MH_US)
C_chisq_MH_AN <- pairwise_c_chisq(MH_AN)
C_chisq_GH_US <- pairwise_c_chisq(GH_US)
C_chisq_GH_AN <- pairwise_c_chisq(GH_AN)
C_chisq_US_AN <- pairwise_c_chisq(US_AN)

C_chisq <- as.data.frame(cbind(C_chisq_MH_GH, C_chisq_MH_US, C_chisq_MH_AN, C_chisq_GH_US, C_chisq_GH_AN, C_chisq_US_AN))

write.table(C_chisq, file="C_chisq_floweringTime.txt", quote=F, sep="\t", row.names=F)


C_chisq_p_MH_GH <- allwise_p_chisq(MH_GH)
C_chisq_p_MH_US <- allwise_p_chisq(MH_US)
C_chisq_p_MH_AN <- allwise_p_chisq(MH_AN)
C_chisq_p_GH_US <- allwise_p_chisq(GH_US)
C_chisq_p_GH_AN <- allwise_p_chisq(GH_AN)
C_chisq_p_US_AN <- allwise_p_chisq(US_AN)

C_chisq_p <- as.data.frame(cbind(C_chisq_p_MH_GH, C_chisq_p_MH_US, C_chisq_p_MH_AN, C_chisq_p_GH_US, C_chisq_p_GH_AN, C_chisq_p_US_AN))

write.table(C_chisq_p, file="C_chisq_p_flowering.txt", quote=F, sep="\t", row.names=F)

#this file contains the candidate gene in flowering time pathway from dong et al. 2012
df <- read.delim("Dong_flowering_outlier_each_pop.txt", header=T)
row.names(df) <- df$geneID
head(df)
df.new <- data.matrix(df[, 2:5])

###null hypotheses: all genes in the pathway have the equal possibility being selected

MH_GH <- df.new[, c(1,3)]
MH_US <- df.new[, 1:2]
MH_AN <- df.new[, c(1,4)]
GH_US <- df.new[, 2:3]
GH_AN <- df.new[, c(3,4)]
US_AN <- df.new[, c(2,4)]

C_chisq_MH_GH <- pairwise_c_chisq(MH_GH)
C_chisq_MH_US <- pairwise_c_chisq(MH_US)
C_chisq_MH_AN <- pairwise_c_chisq(MH_AN)
C_chisq_GH_US <- pairwise_c_chisq(GH_US)
C_chisq_GH_AN <- pairwise_c_chisq(GH_AN)
C_chisq_US_AN <- pairwise_c_chisq(US_AN)

C_chisq <- as.data.frame(cbind(C_chisq_MH_GH, C_chisq_MH_US, C_chisq_MH_AN, C_chisq_GH_US, C_chisq_GH_AN, C_chisq_US_AN))
write.table(C_chisq, file="C_chisq_flowering_dong2012.txt", quote=F, sep="\t", row.names=F)


C_chisq_p_MH_GH <- allwise_p_chisq(MH_GH)
C_chisq_p_MH_US <- allwise_p_chisq(MH_US)
C_chisq_p_MH_AN <- allwise_p_chisq(MH_AN)
C_chisq_p_GH_US <- allwise_p_chisq(GH_US)
C_chisq_p_GH_AN <- allwise_p_chisq(GH_AN)
C_chisq_p_US_AN <- allwise_p_chisq(US_AN)

C_chisq_p <- as.data.frame(cbind(C_chisq_p_MH_GH, C_chisq_p_MH_US, C_chisq_p_MH_AN, C_chisq_p_GH_US, C_chisq_p_GH_AN, C_chisq_p_US_AN))

write.table(C_chisq_p, file="C_chisq_p_flowering_dong2012.txt", quote=F, sep="\t", row.names=F)



