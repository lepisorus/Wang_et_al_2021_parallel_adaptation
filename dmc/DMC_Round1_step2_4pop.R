rec = 0.0000007 #per base pair recombination rate estimate for the region
Ne = 166000
numPops = 4
sampleSizes = c(8, 6, 12, 10)
selPops = c(1, 2, 3, 4)

F_estimate = readRDS("neutralF.RDS")


positions = seq(1, 20001, by=10)

numBins = 1000

selSite = 10001
sels = c(1e-4, 1e-3, 0.01, seq(0.02, 0.1, by = 0.01), seq(0.15, 0.3, by = 0.05), seq(0.4, 0.6, by = 0.1))
times = c(10, 100, 1000, 2000, 4000, 6000, 8000, 10000, 100000)
gs = c(1/(2*Ne), 10^-(4:1))
migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)
sources = selPops

source("genSelMatrices_individualModes.R")

#For model 1 (all selected populations have independent mutations of beneficial allele):
FOmegas_ind = lapply(sels, function(sel) {
  calcFOmegas_indSweeps(sel)
})

saveRDS(FOmegas_ind, "FOmegas_ind_Round1.RDS")

##For model 2 (all selected populations share beneficial allele via migration):
FOmegas_mig = lapply(sels ,function(sel) {
  lapply(migs, function(mig) {
    lapply(sources, function(my.source) {
      calcFOmegas_mig(sel, mig, my.source)
    })
  })
})

saveRDS(FOmegas_mig, "FOmegas_mig_Round1.RDS")

##For model 3  (the beneficial allele was standing in the ancestor of all selected populations):
FOmegas_sv = lapply(sels, function(sel) {
  lapply(gs, function(g) {
    lapply(times, function(time) {    
        calcFOmegas_stdVar(sel, g, time)   
    })
  })
})

saveRDS(FOmegas_sv, "FOmegas_sv_Round1.RDS")


