library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup.RData')
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

time.set <- rep(0, 5)
iters <- 30000
burn <- 25000
update <- 1000
# iters <- 100
# burn <- 50
# update <- 10

method <- "gaussian"
nknots <- 1
threshold <- 0
start <- proc.time()

set.seed(1)

tic.set <- proc.time()
fit.1 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, method=method,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[1] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0
start <- proc.time()

set.seed(2)

tic.set <- proc.time()
fit.2 <- mcmc(y=y, s=s.scale, x=X, method=method,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[2] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0
start <- proc.time()

set.seed(3)

tic.set <- proc.time()
fit.3 <- mcmc(y=y, s=s.scale, x=X, method=method,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[3] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0.90
start <- proc.time()

tic.set <- proc.time()
fit.4 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, method=method,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[4] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0.90
start <- proc.time()

tic.set <- proc.time()
fit.5 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, method=method,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

save(fit.1, fit.2, fit.3, fit.4, fit.5, file="OzoneFull.RData")


