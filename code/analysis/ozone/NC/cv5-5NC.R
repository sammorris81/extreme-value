options(warn=2)
library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup-nc.RData')
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

setting <- 5
method <- "t"
nknots <- 5
threshold <- 0.90
outputfile <- paste("cv5-", setting, "NC.RData", sep="")
start <- proc.time()

fit <- vector(mode="list", length=5)

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
	
	tic.set <- proc.time()
	fit[[val]] <- mcmc(y=y.o, s=S.o, x=X.o, x.pred=X.p, s.pred=S.p,
	                   method=method, skew=F,
	                   thresh=threshold, nknots=nknots, 
                       iters=30000, burn=25000, update=1000, iterplot=F,
                       beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
                       nu.init=0.5, alpha.init=0.5)
	toc.set <- proc.time()
	time.set <- (toc.set - tic.set)[3]
	
	elap.time.val <- (proc.time() - start)[3]
	avg.time.val <- elap.time.val / val
	cat("CV", val, "finished", round(avg.time.val, 2), "per dataset \n")
	save(fit, file=outputfile)
}
