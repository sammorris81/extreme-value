library(fields)
library(geoR)
library(mvtnorm)
library(evd)

load('cv.RData')
source('./densities.R')
source('./mcmc.R')

setting <- 2

threshold <- thresholds[setting]
pred.gpd <- pred.gam <- pred.mvn <- array(NA, dim=c(11, nt, iters, 5))
params <- array(NA, dim=c(iters, 7, 5))
beta.gpd <- beta.gam <- beta.mvn <- array(NA, dim=c(iters, dim(X)[3], 5))

outputfile <- paste("cv-", setting, ".RData", sep="")
start <- proc.time()

for(val in 1:5){
	
	set.seed(setting*100 + val)
	
	cat("CV", val, "started \n")
	val.idx <- cv.lst[[val]]
	y.o <- y[-val.idx,]
	X.o <- X[-val.idx,,]
	S.o <- s.scale[-val.idx,]
	
	y.p <- y[val.idx,]
	X.p <- X[val.idx,,]
	S.p <- s.scale[val.idx,]
	
	tic.gpd <- proc.time()
	fit.gpd <- mcmc(y=y.o, S=S.o, X=X.o, Sp=S.p, Xp=X.p, thresh=threshold,
				    r.model="GPD", iters=iters, burn=burn, update=100,
				    thin=1)
	toc.gpd <- proc.time()
	time.gpd <- (toc.gpd - tic.gpd)[3]
	
	pred.gpd[,,,val] <- fit.gpd$yp
	params[,,val] <- fit.gpd$params
	beta.gpd[,,val] <- fit.gpd$beta
	
	cat("  gpd finished in", time.gpd, "\n")
	
	tic.gam <- proc.time()
	fit.gam <- mcmc(y=y.o, S=S.o, X=X.o, Sp=S.p, Xp=X.p, thresh=threshold,
				    r.model="gamma", iters=iters, burn=burn, update=100,
				    thin=1)
	toc.gam <- proc.time()
	time.gam <- (toc.gam - tic.gam)[3]
	
	pred.gam[,,,val] <- fit.gam$yp
	beta.gam[,,val] <- fit.gam$beta
	
	cat("  gam finished in", time.gam, "\n")
	
	tic.mvn <- proc.time()
	fit.mvn <- mcmc(y=y.o, S=S.o, X=X.o, Sp=S.p, Xp=X.p, thresh=threshold,
				    r.model="fixed", iters=iters, burn=burn, update=100,
				    thin=1)
	toc.mvn <- proc.time()
	time.mvn <- (toc.mvn - tic.mvn)[3]
	
	pred.mvn[,,,val] <- fit.mvn$yp
	beta.mvn[,,val] <- fit.mvn$beta
	
	cat("  fixed finished in", time.mvn, "\n")
	elap.time.val <- (proc.time() - start)[3]
	avg.time.val <- elap.time.val / val
	cat("CV", val, "finished", round(avg.time.val, 2), "per dataset \n")
	save(pred.gpd, pred.gam, pred.mvn, params, beta.gpd, beta.gam, beta.mvn, file=outputfile)
}
