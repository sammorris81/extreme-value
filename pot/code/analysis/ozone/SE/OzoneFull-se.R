library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup-se.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

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
skew <- F
start <- proc.time()

set.seed(1)

tic.set <- proc.time()
fit.1 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, method=method, skew=skew,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[1] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0
skew <- F
start <- proc.time()

set.seed(2)

tic.set <- proc.time()
fit.2 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[2] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0
skew <- F
start <- proc.time()

set.seed(3)

tic.set <- proc.time()
fit.3 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[3] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0.90
skew <- F
start <- proc.time()

tic.set <- proc.time()
fit.4 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[4] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0.90
skew <- F
start <- proc.time()

tic.set <- proc.time()
fit.5 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0
skew <- T
start <- proc.time()

tic.set <- proc.time()
fit.6 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0.9
skew <- T
start <- proc.time()

tic.set <- proc.time()
fit.7 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0
skew <- T
start <- proc.time()

tic.set <- proc.time()
fit.8 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

X <- X[, , c(1, 2, 3)]

method <- "gaussian"
nknots <- 1
threshold <- 0
skew <- F
start <- proc.time()

set.seed(1)

tic.set <- proc.time()
fit.9 <- mcmc(y=y, s=s.scale, x=X, 
            thresh=threshold, nknots=nknots, method=method, skew=skew,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[1] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0
skew <- F
start <- proc.time()

set.seed(2)

tic.set <- proc.time()
fit.10 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[2] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0
skew <- F
start <- proc.time()

set.seed(3)

tic.set <- proc.time()
fit.11 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[3] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0.90
skew <- F
start <- proc.time()

tic.set <- proc.time()
fit.12 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots, 
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[4] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0.90
skew <- F
start <- proc.time()

tic.set <- proc.time()
fit.13 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0
skew <- T
start <- proc.time()

tic.set <- proc.time()
fit.14 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 1
threshold <- 0.9
skew <- T
start <- proc.time()

tic.set <- proc.time()
fit.15 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 5
threshold <- 0
skew <- T
start <- proc.time()

tic.set <- proc.time()
fit.16 <- mcmc(y=y, s=s.scale, x=X, method=method, skew=skew,
            thresh=threshold, nknots=nknots,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]

save(fit.1, fit.2, fit.3, fit.4, fit.5, fit.6, fit.7, fit.8, fit.9, fit.10,
     fit.11, fit.12, fit.13, fit.14, fit.15, fit.16, file="OzoneFull.RData")


