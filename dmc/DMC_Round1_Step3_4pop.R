library("MASS")

## Neutral model
sampleSizes = c(8, 6, 12, 10)
F_estimate = readRDS("neutralF.RDS")

numPops = 4
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 

sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)

det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(det_FOmegas_neutral, "det_FOmegas_neutral_Round1.RDS")

inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(inv_FOmegas_neutral, "inv_FOmegas_neutral_Round1.RDS")



## Model 1
FOmegas_ind = readRDS("FOmegas_ind_Round1.RDS")

det_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
    lapply(sel, function(dist) {
        det(dist)
    })
})
saveRDS(det_FOmegas_ind, "det_FOmegas_ind_Round1.RDS")

inv_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
    lapply(sel, function(dist) {
        ginv(dist)
    })
})
saveRDS(inv_FOmegas_ind, "inv_FOmegas_ind_Round1.RDS")



## Model 2
FOmegas_mig = readRDS("FOmegas_mig_Round1.RDS")

det_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
    lapply(sel, function(mig) {
        lapply(mig, function(source) {
            lapply(source, function(dist) {
                det(dist)
            })
        })
    })
})
saveRDS(det_FOmegas_mig, "det_FOmegas_mig_Round1.RDS")

inv_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
    lapply(sel, function(mig) {
        lapply(mig, function(source) {
            lapply(source, function(dist) {
                ginv(dist)
            })
        })
    })
})
saveRDS(inv_FOmegas_mig, "inv_FOmegas_mig_Round1.RDS")



## Model 3 ##Standing variant without source model
FOmegas_sv = readRDS("FOmegas_sv_Round1.RDS")


det_FOmegas_stdVar = lapply(FOmegas_sv, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(dist) {
				det(dist)
			})
		})
	})
})
saveRDS(det_FOmegas_stdVar, "det_FOmegas_sv_Round1.RDS")

inv_FOmegas_stdVar = lapply(FOmegas_sv, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(dist) {
				ginv(dist)
			})
		})
	})
})
saveRDS(inv_FOmegas_stdVar, "inv_FOmegas_sv_Round1.RDS")

