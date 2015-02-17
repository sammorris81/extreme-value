library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup-se.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

#### Create places for prediciton maps
s1.preds <- seq(1050, 1800, length=25)
s2.preds <- seq(-860, -250, length=25)
s.preds <- expand.grid(s1.preds, s2.preds)
s.preds <- s.preds[(s.preds[, 2] >= (1.33 * s.preds[, 1] - 2815)), ]  # atlantic
s.preds <- s.preds[(s.preds[, 2] < (0.75 * s.preds[, 1] - 1285)), ]  # tennessee
s.preds <- s.preds[(s.preds[, 2] >= (-6.2 * s.preds[, 1] + 5960)), ]  # tennessee

# plot(s.preds)
# lines(l)

s.scale.preds <- matrix(NA, nrow=nrow(s.preds), ncol=ncol(s.preds))
s.scale.preds[,1] <- (s.preds[,1] - range(s[,1])[1])/(range(s[,1])[2] - range(s[,1])[1])
s.scale.preds[,2] <- (s.preds[,2] - range(s[,2])[1])/(range(s[,2])[2] - range(s[,2])[1])

X <- array(1, dim=c(ns, nt, 6))
X.preds <- array(1, dim=c(nrow(s.preds), nt, 6))
for (t in 1:nt) {
  X[, t, 2] <- s.scale[, 1]
  X[, t, 3] <- s.scale[, 2]
  X[, t, 4] <- s.scale[, 1]^2
  X[, t, 5] <- s.scale[, 2]^2
  X[, t, 6] <- s.scale[, 1] * s.scale[, 2]
  
  X.preds[, t, 2] <- s.scale.preds[, 1]
  X.preds[, t, 3] <- s.scale.preds[, 2]
  X.preds[, t, 4] <- s.scale.preds[, 1]^2
  X.preds[, t, 5] <- s.scale.preds[, 2]^2
  X.preds[, t, 6] <- s.scale.preds[, 1] * s.scale.preds[, 2]
}

time.set <- rep(0, 5)
iters <- 30000
burn <- 25000
update <- 500


method <- "t"
nknots <- 1
threshold <- 0.90
skew <- F
start <- proc.time()

set.seed(2)

tic.set <- proc.time()
fit.2 <- mcmc(y=y, s=s.scale, x=X, x.pred=X.preds, s.pred=s.scale.preds,
            thresh=threshold, nknots=nknots, method=method, skew=skew,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[2] <- (toc.set - tic.set)[3]


save(fit.2, file="OzoneFull2.RData")


