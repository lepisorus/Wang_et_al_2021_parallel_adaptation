#! All new comments are marked with "#!"


#Calculate composite likelihoods


#frequency file
df <- read.table("MH.GH.US.AN.REFfrq", header=F) #! this will change once only have file with the 4 pops and invariable sites removed

head(df)

##setting up parameters
numPops = 4 #! change this to 4
#these values must be same as used to calculate Calculate F^(S) matrices above
numBins = 1000

rec = 0.0000007 #per base pair recombination rate estimate for the region
Ne = 166000
sels = c(1e-4, 1e-3, 0.01, seq(0.02, 0.1, by = 0.01), seq(0.15, 0.3, by = 0.05), seq(0.4, 0.6, by = 0.1)) #! reduce this for speed up
times = c(10, 100, 1000, 2000, 4000, 6000, 8000, 10000, 100000)
gs = c(1/(2*Ne), 10^-(4:1))
migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)
selPops = c(1, 2, 3, 4) #! change this to match indices of AH and MH
sources = selPops #! this will mean that either AH or MH is able to be the source of the beneficial allele in the migration model


#! Randomize allele frequencies wrt reference allele and mean-center here (outside of the loop)
freqs_notRand <- t(df[, 3:6]) #! will want different columns once this file only contains 4 populations

randFreqs = apply(freqs_notRand, 2, function(my.freqs) {
  if(runif(1) < 0.5) {
    my.freqs = 1 - my.freqs
  }
  my.freqs
})

saveRDS(randFreqs, "round1.randFreqs") 
#! save these in case we need to re-run analyses (then skip previous 2 lines and simply re-load data)
#! want to use same randomization throughout (i.e. only randomize once and use same frequencies for all analyses)

freqs_all = t(randFreqs)

epsilons_all = rowMeans(freqs_all)

M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 
freqs_MC_all = sapply(1 : nrow(freqs_all), function(i) Tmatrix %*% freqs_all[i,])
#freqs_MC_all <- matrix(freqs_MC_all, nrow=M-1)

## how to loop through positions, freqs and selSites

focalRegions <- read.delim("fourHighOutlierSNP.focalRegions.txt", header=F)
names(focalRegions) <- c("chr", "start", "end")

likelihood <- data.frame(chr=as.numeric(), selSite=as.numeric(), neutral=as.numeric(), ind=as.numeric(), mig=as.numeric(), sv=as.numeric(), nSNP=as.numeric(), time10=as.numeric(), time100=as.numeric(), time1000=as.numeric(), time2000=as.numeric(), time4000=as.numeric(), time6000=as.numeric(), time8000=as.numeric(), time10000=as.numeric(), time100000=as.numeric(), migSource=as.numeric(), selCoeff_ind=as.numeric(), selCoeff_mig=as.numeric(), selCoeff_sv=as.numeric())

for (row in 1:nrow(focalRegions)) {
    chr <- focalRegions[row, "chr"]
    start <- focalRegions[row, "start"]
    end <- focalRegions[row, "end"]
    selSite = start + 10000

##one could test the script by running through the first focal region, without the for loop

    #chr <- 10
    #start <- 743882
    #end <- 763882
    #selSite = 753882    
    
    
    #remove selected site
    which.rows = which(df$V1 == chr & df$V2 >= start & df$V2 <= end & df$V2 != selSite)
    positions = df$V2[which.rows]
    nSNP = length(positions)
  if (nSNP >= 2) {  
    #freq
    freqs_MC = freqs_MC_all[ , which.rows]
    epsilons = epsilons_all[which.rows]
#    freqs_MC = matrix(freqs_MC, nrow=M-1)
   
   source("calcCompositeLike_1site.R") #! changed
    #calculate composite likelihoods
    
 ## Neutral model
det_FOmegas_neutral = readRDS("det_FOmegas_neutral_Round1.RDS")
inv_FOmegas_neutral = readRDS("inv_FOmegas_neutral_Round1.RDS")
compLikelihood_neutral = calcCompLikelihood_neutral.1site(det_FOmegas_neutral, inv_FOmegas_neutral) #! changed
saveRDS(compLikelihood_neutral, paste("compLikelihood_neutral_Round1_", chr, "_", selSite, ".RDS", sep=""))

## Model 1
det_FOmegas_ind = readRDS("det_FOmegas_ind_Round1.RDS")
inv_FOmegas_ind = readRDS("inv_FOmegas_ind_Round1.RDS")
compLikelihood_ind = lapply(1 : length(sels), function(sel) calcCompLikelihood_1par.1site(det_FOmegas_ind, inv_FOmegas_ind, sel)) #! changed
saveRDS(compLikelihood_ind, paste("compLikelihood_ind_Round1_", chr, "_", selSite, ".RDS", sep=""))
mcleS_ind = sels[which.max(unlist(lapply(1 : length(sels), function(i) max(unlist(compLikelihood_ind[[i]])))))]

## Model 2
det_FOmegas_mig = readRDS("det_FOmegas_mig_Round1.RDS")
inv_FOmegas_mig = readRDS("inv_FOmegas_mig_Round1.RDS")
compLikelihood_mig = lapply(1 : length(sels), function(sel) {
        lapply(1 : length(migs), function(mig) {
            lapply(1 : length(sources), function(my.source) {
                calcCompLikelihood_3par.1site(det_FOmegas_mig, inv_FOmegas_mig, sel, mig,
                                        my.source)
        })
    })
}) #! changed
saveRDS(compLikelihood_mig, paste("compLikelihood_mig_Round1_", chr, "_", selSite, ".RDS", sep=""))
max_source1 = max(unlist(lapply(1: length(compLikelihood_mig), function(i) unlist(compLikelihood_mig[[i]][[1]]))))
max_source2 = max(unlist(lapply(1: length(compLikelihood_mig), function(i) unlist(compLikelihood_mig[[i]][[2]]))))
max_source3 = max(unlist(lapply(1: length(compLikelihood_mig), function(i) unlist(compLikelihood_mig[[i]][[3]]))))
max_source4 = max(unlist(lapply(1: length(compLikelihood_mig), function(i) unlist(compLikelihood_mig[[i]][[4]]))))
mcleSource = which.max(c(max_source1, max_source2，max_source3，max_source4))
mcleS_mig = sels[which.max(unlist(lapply(1 : length(sels), function(i) max(unlist(compLikelihood_mig[[i]])))))]

## Model 3
det_FOmegas_sv = readRDS("det_FOmegas_sv_Round1.RDS")
inv_FOmegas_sv = readRDS("inv_FOmegas_sv_Round1.RDS")
compLikelihood_sv = lapply(1 : length(sels), function(sel) {
        lapply(1 : length(gs), function(g) {
            lapply(1 : length(times), function(t) {
                calcCompLikelihood_3par.1site(det_FOmegas_sv, inv_FOmegas_sv, sel, g, t)
        })
    })
}) #! changed
saveRDS(compLikelihood_sv, paste("compLikelihood_sv_Round1_", chr, "_", selSite, ".RDS", sep=""))
compLikelihood_sv_tFirst = lapply(1 : length(times), function(time) sapply(1: length(sels), function(i) sapply(1 : length(gs), function(j) compLikelihood_sv[[i]][[j]][[time]])))
mcleT = sapply(1: length(compLikelihood_sv_tFirst), function(i) max(unlist(compLikelihood_sv_tFirst[[i]])))
mcleS_sv = sels[which.max(unlist(lapply(1 : length(sels), function(i) max(unlist(compLikelihood_sv[[i]])))))]


#likelihood[nrow(likelihood) + 1,] = list(chr=chr, selSite=selSite, neutral=unlist(composite_likelihood_neutral), ind=unlist(composite_likelihood_ind), mig=unlist(composite_likelihood_mig), sv=unlist(composite_likelihood_sv))
 
likelihood[nrow(likelihood) + 1,] = list(chr=chr, selSite=selSite, neutral=max(unlist(compLikelihood_neutral)), ind=max(unlist(compLikelihood_ind)), mig=max(unlist(compLikelihood_mig)), sv=max(unlist(compLikelihood_sv)), nSNP=nSNP, time10=mcleT[1], time100=mcleT[2], time1000=mcleT[3], time2000=mcleT[4], time4000=mcleT[5], time6000=mcleT[6], time8000=mcleT[7], time10000=mcleT[8], time100000=mcleT[9], migSource=mcleSource, selCoeff_ind=mcleS_ind, selCoeff_mig=mcleS_mig, selCoeff_sv=mcleS_sv)
#! I commented out the line you wrote and modified slightly because I figured you just wanted the maximum composite likelihood across all parameters to summarize

}#if
 }#for

write.table(likelihood, file="likelihoodSummary.txt", sep="\t", row.names=F, quote=F)

