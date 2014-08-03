library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup-se.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

#### Create places for prediciton maps
s1.preds <- seq(1050, 1800, length=50)
s2.preds <- seq(-860, -250, length=50)
s.preds <- expand.grid(s1.preds, s2.preds)
s.preds <- s.preds[(s.preds[, 2] >= (1.33 * s.preds[, 1] - 2815)), ]  # atlantic
s.preds <- s.preds[(s.preds[, 2] < (0.75 * s.preds[, 1] - 1285)), ]  # tennessee
s.preds <- s.preds[(s.preds[, 2] >= (-6.2 * s.preds[, 1] + 5960)), ]  # tennessee

# plot(s.preds)
# lines(l)

s.scale.preds <- matrix(NA, nrow=nrow(s.preds), ncol=ncol(s.preds))
s.scale.preds[,1] <- (s.preds[,1] - range(s[,1])[1])/(range(s[,1])[2] - range(s[,1])[1])
s.scale.preds[,2] <- (s.preds[,2] - range(s[,2])[1])/(range(s[,2])[2] - range(s[,2])[1])

x.preds[, ] <- cbind(rep(1, nrow(s.preds)), s.scale.preds)
X <- X[, , c(1, 2, 3)]
X.preds <- array(1, dim=c(nrow(s.preds), nt, 3))
for (t in 1:nt) {
  X.preds[, t, 2] <- s.scale.preds[, 1]
  X.preds[, t, 3] <- s.scale.preds[, 2]
}


time.set <- rep(0, 5)
iters <- 30000
burn <- 25000
update <- 1000

method <- "gaussian"
nknots <- 1
threshold <- 0
skew <- F
start <- proc.time()

set.seed(1)

tic.set <- proc.time()
fit.1 <- mcmc(y=y, s=s.scale, x=X, x.pred=X.preds, s.pred=s.scale.preds,
            thresh=threshold, nknots=nknots, method=method, skew=skew,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[1] <- (toc.set - tic.set)[3]

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

method <- "t"
nknots <- 1
threshold <- 0.90
skew <- T
start <- proc.time()

set.seed(3)

tic.set <- proc.time()
fit.3 <- mcmc(y=y, s=s.scale, x=X, x.pred=X.preds, s.pred=s.scale.preds,
            thresh=threshold, nknots=nknots, method=method, skew=skew,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[3] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 3
threshold <- 0.90
skew <- F
start <- proc.time()

set.seed(4)

tic.set <- proc.time()
fit.4 <- mcmc(y=y, s=s.scale, x=X, x.pred=X.preds, s.pred=s.scale.preds,
            thresh=threshold, nknots=nknots, method=method, skew=skew,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[4] <- (toc.set - tic.set)[3]

method <- "t"
nknots <- 3
threshold <- 0.90
skew <- T
start <- proc.time()

set.seed(5)

tic.set <- proc.time()
fit.5 <- mcmc(y=y, s=s.scale, x=X, x.pred=X.preds, s.pred=s.scale.preds,
            thresh=threshold, nknots=nknots, method=method, skew=skew,
            iters=iters, burn=burn, update=update, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set[5] <- (toc.set - tic.set)[3]


save(fit.1, fit.2, fit.3, fit.4, fit.5, file="OzoneFull.RData")

# make the plots
library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

#### Setup from Brian
load("../OzoneData.RData")
y <- CMAQ_OZONE$y
x <- CMAQ_OZONE$x
s <- CMAQ_OZONE$s
l <- CMAQ_OZONE$poly

#Extract subset
# NC<-(s[,1]>xb[1]) & (s[,1]<xb[2]) & (s[,2]>yb[1]) & (s[,2]<yb[2])
AL <- 1:26
DE <- 27:30
DC <- 31:33
FL <- 34:45
GA <- 46:69
IL <- 70:76
IN <- 77:104
KY <- 105:135
MD <- 136:150
MS <- 151:154
NJ <- 155
NC <- 156:196
OH <- 197:223
PA <- 224:231
SC <- 232:252
TN <- 253:274
VA <- 275:299
WV <- 300:307

SE <- c(GA, NC, SC)
s<-s[SE,]
x<-x[SE,]
y<-y[SE,]
plot(s)
lines(l)

load('./OzoneFull2.RData')
yp <- fit.2$yp
q.95 <- apply(yp, 2, quantile, probs=c(0.95))

p.exceed.75 <- rep(NA, 439)
# for each site, how much of the posterior predictive is above 75ppb
for (i in 1:439) {
  p.exceed.75[i] <- mean(yp[, i, ] > 75)
}