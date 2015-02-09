options(warn=2)
library(fields)
library(SpatialTools)
library(mvtnorm)

rm(list=ls())
load('us-all-setup.RData')
source('../../../R/mcmc.R', chdir=T)
source('../../../R/auxfunctions.R')
source('../max-stab/gen_data.R')
source('../max-stab/MCMC4MaxStable.R')

setting <- 26
outputfile <- paste("us-all-", setting, ".RData", sep="")

start <- proc.time()
thresh <- 75
fit <- vector(mode="list", length=2)

for(val in 1:2){

	set.seed(setting * 100 + val)

	cat("CV", val, "started \n")
	val.idx <- cv.lst[[val]]
	y.o <- t(Y[-val.idx,])
	X.o <- t(X[-val.idx, ,2])
	S.o <- S[-val.idx,]
	knots <- S.o

	y.p <- t(Y[val.idx,])
	X.p <- t(X[val.idx, , 2])
	S.p <- S[val.idx,]

	tic.set <- proc.time()
	fit[[val]] <- maxstable(y=y.o, x=X.o, s=S.o, sp=S.p, xp=X.p, thresh=thresh,
	                        knots=knots, iters=30000, burn=25000, update=500, thin=1)
	toc.set <- proc.time()
	time.set <- (toc.set - tic.set)[3]

	elap.time.val <- (proc.time() - start)[3]
	avg.time.val <- elap.time.val / val
	cat("CV", val, "finished", round(avg.time.val, 2), "per dataset \n")
	save(fit, file=outputfile)
}
