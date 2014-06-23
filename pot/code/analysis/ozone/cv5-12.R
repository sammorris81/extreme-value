library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup.RData')
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

setting <- 12
threshold <- 0.99
nknots <- 5
outputfile <- paste("cv-", setting, ".RData", sep="")
start <- proc.time()

fit <- vector(mode="list", length=5)

for(val in 1:5){
	
	set.seed(setting*100 + val)
	
	cat("CV", val, "started \n")
	val.idx <- cv.lst[[val]]
	y.o <- y[-val.idx,]
	X.o <- X[-val.idx,,]
	S.o <- s[-val.idx,]
	
	rho.ml <- RhoML(S.o, y.o)
	
	y.p <- y[val.idx,]
	X.p <- X[val.idx,,]
	S.p <- s[val.idx,]
	
	tic.set <- proc.time()
	fit[[val]] <- mcmc(y=y.o, s=S.o, x=X.o, x.pred=X.p, s.pred=S.p, 
	                   thresh=threshold, nknots=nknots, 
                       iters=30000, burn=25000, update=1000, iterplot=F,
                       beta.init=beta.init, tau.init=tau.init, 
                       rho.init=rho.ml, fixrho=T,
                       nu.init=0.5, alpha.init=0.5, delta.init=0, fixdelta=T,
                       fixz=T, z.init=matrix(0, nrow=nknots, ncol=nt), scale=F)
	toc.set <- proc.time()
	time.set <- (toc.set - tic.set)[3]
	
	elap.time.val <- (proc.time() - start)[3]
	avg.time.val <- elap.time.val / val
	cat("CV", val, "finished", round(avg.time.val, 2), "per dataset \n")
	save(fit, file=outputfile)
}
