library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup.RData')
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

threshold <- 0
nknots <- 1
start <- proc.time()

set.seed(1)

tic.set <- proc.time()
fit.1 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, 
            iters=30000, burn=25000, update=1000, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5, delta.init=0, fixdelta=T, scale=T,
            fixz=T, z.init=matrix(rep(0, nt), 1, nt))
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

threshold <- 0
nknots <- 5
start <- proc.time()

set.seed(2)

tic.set <- proc.time()
fit.5 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, 
            iters=30000, burn=25000, update=1000, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5, delta.init=0, fixdelta=T, scale=T,
            fixz=T, z.init=matrix(rep(0, nt * nknots), nknots, nt))
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

threshold <- 0
nknots <- 10
start <- proc.time()

set.seed(3)

tic.set <- proc.time()
fit.10 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, 
            iters=30000, burn=25000, update=1000, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5, delta.init=0, fixdelta=T, scale=T,
            fixz=T, z.init=matrix(rep(0, nt * nknots), nknots, nt))
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

save(fit.1, fit.5, fit.10, file="OzoneFull.RData")


