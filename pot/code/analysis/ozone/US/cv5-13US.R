options(warn=2)
library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup-us.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

setting <- 13
method <- "t"
nknots <- 10
keep.knots <- F
threshold <- 90
thresh.quant <- F
skew <- T
outputfile <- paste("cv5-", setting, "US.RData", sep="")

start <- proc.time()

fit <- vector(mode="list", length=5)

for(val in 1:5){
	
	set.seed(setting*100 + val)
	
	cat("CV", val, "started \n")
	val.idx <- cv.lst[[val]]
	y.o <- aqs[-val.idx,]
	X.o <- X[-val.idx,,]
	S.o <- s[-val.idx,]
	
	y.p <- aqs[val.idx,]
	X.p <- X[val.idx,,]
	S.p <- s[val.idx,]
	
	tic.set <- proc.time()
	fit[[val]] <- mcmc(y=y.o, s=S.o, x=X.o, x.pred=X.p, s.pred=S.p,
	                   method=method, skew=skew, keep.knots=keep.knots,
	                   thresh.all=threshold, thresh.quant=thresh.quant, nknots=nknots, 
                       iters=30000, burn=25000, update=500, iterplot=F,
                       beta.init=beta.init, tau.init=tau.init, rho.init=1000,
                       nu.init=0.5, alpha.init=0.5)
	toc.set <- proc.time()
	time.set <- (toc.set - tic.set)[3]
	
	elap.time.val <- (proc.time() - start)[3]
	avg.time.val <- elap.time.val / val
	cat("CV", val, "finished", round(avg.time.val, 2), "per dataset \n")
	save(fit, file=outputfile)
}
