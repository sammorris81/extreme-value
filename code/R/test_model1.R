# clean out previous variables
rm(list=ls())

# libraries
library(fields)
library(SpatialTools)
library(emulator)

# necessary functions
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')

# data settings
beta.t      <- c(10, 0, 0)
nu.t        <- 0.5
gamma.t     <- 0.9
rho.t       <- 1
tau.alpha.t <- 3
tau.beta.t  <- 8
lambda.t    <- 2
lambda.1.t  <- sign(lambda.t)  # direction of skew
lambda.2.t  <- 1 / lambda.t^2

# covariate data
set.seed(20)
s         <- cbind(runif(144, 0, 10), runif(144, 0, 10))
knots.x   <- seq(1, 9, length=12)
knots.gev <- expand.grid(knots.x, knots.x)
ns        <- nrow(s)
nt        <- 50
nsets     <- 50
nsettings <- 7
ntest     <- 44

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

d <- rdist(s)
diag(d) <- 0
cor.t <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho.t),
                       error.var=0, smoothness=nu.t, finescale.var=0)
C <- gamma.t * cor.t
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec.t <- CC$prec
logdet.prec.t <- CC$logdet.prec

################################################################################
# MCMC Tests
################################################################################

# Test 1a: No-skew, 1 knot, no time series
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=0, phi.z=0, phi.w=0, phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=10000, burn=5000, update=100,
            thresh.all=0, skew=FALSE, nknots=1, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
# RESULTS: PASS

# Test 1b: No-skew, 1 knot, no time series
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=0, phi.z=0, phi.w=0, phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=10000, burn=5000, update=100, lambda.init=1,
            thresh.all=0, skew=TRUE, nknots=1, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

# Test 2: Skew, 1 knot, no time series
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(20)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=2, phi.z=0, phi.w=0, phi.tau=0)
# z1: 3.764, z10=3.784, z21=3.67

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=10000, burn=5000, update=100,
            thresh.all=0, skew=TRUE, nknots=1, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
# RESULTS: PASS

# Test 3: Non-Skew, 3 knots, no time series
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(20)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=3, lambda=0, phi.z=0, phi.w=0, phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100,
            thresh.all=0, skew=FALSE, nknots=3, rho.upper=15, nu.upper=10,
            # fixknots=TRUE, knots.init=data$knots,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
# RESULTS: PASS

# Test 4a: Skew, 3 knots, no time series
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS1(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=3, lambda=-2, phi.z=0, phi.w=0, phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100, tau.init=0.375,
            thresh.all=0, skew=TRUE, nknots=3, rho.upper=15, nu.upper=10,
            # knots.init=data$knots, fixknots=TRUE,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
# RESULTS: PASS (as long you specify a reasonable bounding box for s)
# Although it eventually starts to get near the right values for most of the
# parameters, it can take upwards of 5000 - 6000 iterations
# Another minor tweak to the MCMC, and it would appear that starting the knots
# in the center of the space tends to stabilize a bit more quickly (2000-3000
# iterations)
# Can't tell for sure, but it looks like it may be a slight improvement to
# adjust the MH standard deviation a little less often

# Test 4b: Skew, 3 knots, no time series, thresholding at 0.80
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS1(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=3, lambda=2, phi.z=0, phi.w=0, phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100, tau.init=0.375,
            thresh.all=0.80, skew=FALSE, nknots=3, rho.upper=15, nu.upper=10,
            # knots.init=data$knots, fixknots=TRUE,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

# Test 5: Non-skew, 1 knots, time series on tau
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=0, phi.z=0, phi.w=0, phi.tau=0.8)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100,
            thresh.all=0, skew=FALSE, nknots=1, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=TRUE, temporalz=FALSE)
# RESULTS: PASS

# Test 6: Skew, 1 knots, time series on z, no time series on tau
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=2, phi.z=0.8, phi.w=0, phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100, lambda.init=0,
            thresh.all=0, skew=TRUE, nknots=1, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=FALSE, temporalz=TRUE)
# RESULTS: PASS

# Test 7: Skew, 1 knots, time series on z and tau
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=2, phi.z=0.8, phi.w=0,
                   phi.tau=0.5)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100,
            thresh.all=0, skew=TRUE, nknots=1, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=FALSE, temporaltau=TRUE, temporalz=TRUE)
# RESULTS: PASS

# Test 8: Non-skew, 3 knots, time series on tau and knots
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=3, lambda=0, phi.z=0, phi.w=0.9,
                   phi.tau=0.8)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100,
            thresh.all=0, skew=FALSE, nknots=3, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=TRUE, temporaltau=TRUE, temporalz=FALSE)
# RESULTS: PASS

# Test 9: Skew, 3 knots, time series on z and knots, no time series on tau
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=3, lambda=2, phi.z=0.8, phi.w=0.9,
                   phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100,
            thresh.all=0, skew=TRUE, nknots=3, rho.upper=15, nu.upper=10,
            # fixknots=TRUE, knots.init=data$knots,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=TRUE, temporaltau=FALSE, temporalz=TRUE)
# RESULTS: PASS (as long as knots start out close to the middle).
# Note: The time series parameters don't appear to get exactly to the correct
# spot, but the rest of the parameters do reasonably well.
# I think this may be related to lack of identifiability within a partition.
# Rho, nu, gamma and lambda also are in the correct neighborhood,
# but not exactly right

# Test 10: Skew, 3 knots, time series on z, knots, and tau
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=3, lambda=2, phi.z=0.8, phi.w=0.9,
                   phi.tau=0.8)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE,
            iterplot=T, iters=7000, burn=5000, update=100,
            thresh.all=0, skew=TRUE, nknots=3, rho.upper=15, nu.upper=10,
            min.s=c(0, 0), max.s=c(10, 10),
            temporalw=TRUE, temporaltau=TRUE, temporalz=TRUE)
# RESULTS: PASS

# Test 11: Predictions
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')

# storage for predictions
fit.1 <- fit.2 <- data <- vector("list", length=5)
for (i in 1:5) {
  set.seed(i)
  data[[i]] <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                          rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                          dist="t", nknots=3, lambda=3, phi.z=0.8, phi.w=0.8,
                          phi.tau=0.8)

  s.o <- s[1:100, ]
  x.o <- x[1:100, , ]
  y.o <- data[[i]]$y[1:100, ]
  s.p <- s[101:144, ]
  x.p <- x[101:144, , ]
  y.p <- data[[i]]$y[101:144, ]

  fit.1[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p, method="t",
                     thresh.quant=TRUE, iterplot=T, iters=10000, burn=7000,
                     update=100, thresh.all=0, skew=TRUE, nknots=3,
                     rho.upper=15, nu.upper=10, min.s=c(0, 0), max.s=c(10, 10),
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  fit.2[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p, method="t",
                     thresh.quant=TRUE, iterplot=T, iters=10000, burn=7000,
                     update=100, thresh.all=0, skew=TRUE, nknots=3,
                     rho.upper=15, nu.upper=10, min.s=c(0, 0), max.s=c(10, 10),
                     temporalw=TRUE, temporaltau=TRUE, temporalz=TRUE)
}

# come up with quantile and brier scores
probs <- seq(0.90, 0.99, by=0.01)
quant.scores.ts  <- brier.scores.ts  <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.nts <- brier.scores.nts <- matrix(NA, nrow=5, ncol=length(probs))
for (i in 1:5) {
  threshs <- quantile(data[[i]]$y, probs=probs)
  quant.scores.ts[i, ]  <- QuantScore(preds=fit.1[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.nts[i, ] <- QuantScore(preds=fit.2[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.ts[i, ]  <- BrierScore(preds=fit.1[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.nts[i, ] <- BrierScore(preds=fit.2[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
}

# one bs and qs for each method
quant.score <- brier.score <- matrix(NA, nrow=2, ncol=length(probs))
quant.score[1, ] <- apply(quant.scores.ts, 2, mean)
quant.score[2, ] <- apply(quant.scores.nts, 2, mean)
brier.score[1, ] <- apply(brier.scores.ts, 2, mean)
brier.score[2, ] <- apply(brier.scores.nts, 2, mean)

# Test 12: Skew, 3 knots, time series on z, knots, and tau with thresholding
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=3, lambda=2, phi.z=0, phi.w=0,
                   phi.tau=0)

fit <- mcmc(y=data$y, s=s, x=x, method="t", thresh.quant=TRUE, iterplot=TRUE,
            iters=1000, burn=900, update=100, thresh.all=0.80, skew=FALSE,
            nknots=3, rho.upper=15, nu.upper=10, min.s=c(0, 0), max.s=c(10, 10),
            nknots=3, temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
# RESULTS:

# Test 13: Predictions - 1 knot, skew
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')

# storage for predictions
fit.1 <- fit.2 <- fit.3 <- data <- vector("list", length=5)
for (i in 1:5) {
  set.seed(i)
  data[[i]] <- rpotspatTS1(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                          rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                          dist="t", nknots=1, lambda=3, phi.z=0, phi.w=0,
                          phi.tau=0)

  s.o <- s[1:100, ]
  x.o <- x[1:100, , ]
  y.o <- data[[i]]$y[1:100, ]
  s.p <- s[101:144, ]
  x.p <- x[101:144, , ]
  y.p <- data[[i]]$y[101:144, ]

  cat("Set", i, "Test 13 - fit.1 \n")
  fit.1[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="gaussian", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     skew=FALSE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 13 - fit.2 \n")
  fit.2[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10, lambda.init=0,
                     tau.init=0.375,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 13 - fit.3 \n")
  fit.3[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     #fixknots=TRUE, knots.init=data[[i]]$knots,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=5,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
}

# come up with quantile and brier scores
probs <- seq(0.90, 0.99, by=0.01)
quant.scores.gau <- brier.scores.gau <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.1st <- brier.scores.1st <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.3st <- brier.scores.3st <- matrix(NA, nrow=5, ncol=length(probs))
for (i in 1:5) {
  threshs <- quantile(data[[i]]$y, probs=probs)
  quant.scores.gau[i, ] <- QuantScore(preds=fit.1[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.1st[i, ] <- QuantScore(preds=fit.2[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.3st[i, ] <- QuantScore(preds=fit.3[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.gau[i, ] <- BrierScore(preds=fit.1[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.1st[i, ] <- BrierScore(preds=fit.2[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.3st[i, ] <- BrierScore(preds=fit.3[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
}

# one bs and qs for each method
quant.score <- brier.score <- matrix(NA, nrow=3, ncol=length(probs))
quant.score[1, ] <- apply(quant.scores.gau, 2, mean)
quant.score[2, ] <- apply(quant.scores.1st, 2, mean)
quant.score[3, ] <- apply(quant.scores.3st, 2, mean)
brier.score[1, ] <- apply(brier.scores.gau, 2, mean)
brier.score[2, ] <- apply(brier.scores.1st, 2, mean)
brier.score[3, ] <- apply(brier.scores.3st, 2, mean)

# Test 14: Predictions - 5 knots, skew
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')

# storage for predictions
fit.1 <- fit.2 <- fit.3 <- data <- vector("list", length=5)
for (i in 1:5) {
  set.seed(i + 5)
  data[[i]] <- rpotspatTS1(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                          rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                          dist="t", nknots=5, lambda=3, phi.z=0, phi.w=0,
                          phi.tau=0)

  s.o <- s[1:100, ]
  x.o <- x[1:100, , ]
  y.o <- data[[i]]$y[1:100, ]
  s.p <- s[101:144, ]
  x.p <- x[101:144, , ]
  y.p <- data[[i]]$y[101:144, ]

  cat("Set", i, "Test 14 - fit.1 \n")
  fit.1[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="gaussian", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     skew=FALSE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 14 - fit.2 \n")
  fit.2[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 14 - fit.1 \n")
  fit.3[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     #fixknots=TRUE, knots.init=data[[i]]$knots,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=5,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
}

# come up with quantile and brier scores
probs <- seq(0.90, 0.99, by=0.01)
quant.scores.gau <- brier.scores.gau <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.1st <- brier.scores.1st <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.3st <- brier.scores.3st <- matrix(NA, nrow=5, ncol=length(probs))
for (i in 1:5) {
  threshs <- quantile(data[[i]]$y, probs=probs)
  quant.scores.gau[i, ] <- QuantScore(preds=fit.1[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.1st[i, ] <- QuantScore(preds=fit.2[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.3st[i, ] <- QuantScore(preds=fit.3[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.gau[i, ] <- BrierScore(preds=fit.1[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.1st[i, ] <- BrierScore(preds=fit.2[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.3st[i, ] <- BrierScore(preds=fit.3[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
}

# one bs and qs for each method
quant.score <- brier.score <- matrix(NA, nrow=3, ncol=length(probs))
quant.score[1, ] <- apply(quant.scores.gau, 2, mean)
quant.score[2, ] <- apply(quant.scores.1st, 2, mean)
quant.score[3, ] <- apply(quant.scores.3st, 2, mean)
brier.score[1, ] <- apply(brier.scores.gau, 2, mean)
brier.score[2, ] <- apply(brier.scores.1st, 2, mean)
brier.score[3, ] <- apply(brier.scores.3st, 2, mean)

# Test 15: Predictions - 1 knot, gaussian
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')

# storage for predictions
fit.1 <- fit.2 <- fit.3 <- data <- vector("list", length=5)
for (i in 1:5) {
  set.seed(i + 10)
  data[[i]] <- rpotspatTS1(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                          rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                          dist="gaussian", nknots=1, lambda=0, phi.z=0, phi.w=0,
                          phi.tau=0)

  s.o <- s[1:100, ]
  x.o <- x[1:100, , ]
  y.o <- data[[i]]$y[1:100, ]
  s.p <- s[101:144, ]
  x.p <- x[101:144, , ]
  y.p <- data[[i]]$y[101:144, ]

  cat("Set", i, "Test 15 - fit.1 \n")
  fit.1[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="gaussian", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     skew=FALSE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 15 - fit.2 \n")
  fit.2[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10, lambda.init=0.01,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 15 - fit.3 \n")
  fit.3[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10, lambda.init=0.01,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=5,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
}

# come up with quantile and brier scores
probs <- seq(0.90, 0.99, by=0.01)
quant.scores.gau <- brier.scores.gau <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.1st <- brier.scores.1st <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.3st <- brier.scores.3st <- matrix(NA, nrow=5, ncol=length(probs))
for (i in 1:5) {
  threshs <- quantile(data[[i]]$y, probs=probs)
  quant.scores.gau[i, ] <- QuantScore(preds=fit.1[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.1st[i, ] <- QuantScore(preds=fit.2[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.3st[i, ] <- QuantScore(preds=fit.3[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.gau[i, ] <- BrierScore(preds=fit.1[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.1st[i, ] <- BrierScore(preds=fit.2[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.3st[i, ] <- BrierScore(preds=fit.3[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
}

# one bs and qs for each method
quant.score <- brier.score <- matrix(NA, nrow=3, ncol=length(probs))
quant.score[1, ] <- apply(quant.scores.gau, 2, mean)
quant.score[2, ] <- apply(quant.scores.1st, 2, mean)
quant.score[3, ] <- apply(quant.scores.3st, 2, mean)
brier.score[1, ] <- apply(brier.scores.gau, 2, mean)
brier.score[2, ] <- apply(brier.scores.1st, 2, mean)
brier.score[3, ] <- apply(brier.scores.3st, 2, mean)

# Test 16: Predictions - 1 knot, gaussian
source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')

# storage for predictions
fit.1 <- fit.2 <- fit.3 <- data <- vector("list", length=5)
for (i in 1:5) {
  set.seed(i + 15)
  data[[i]] <- rpotspatTS1(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                          rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                          dist="t", nknots=5, lambda=0, phi.z=0, phi.w=0,
                          phi.tau=0)

  s.o <- s[1:100, ]
  x.o <- x[1:100, , ]
  y.o <- data[[i]]$y[1:100, ]
  s.p <- s[101:144, ]
  x.p <- x[101:144, , ]
  y.p <- data[[i]]$y[101:144, ]

  cat("Set", i, "Test 16 - fit.1 \n")
  fit.1[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="gaussian", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     skew=FALSE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 16 - fit.2 \n")
  fit.2[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=1,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)

  cat("Set", i, "Test 16 - fit.3 \n")
  fit.3[[i]] <- mcmc(y=y.o, s=s.o, x=x.o, x.pred=x.p, s.pred=s.p,
                     method="t", thresh.quant=TRUE, iterplot=T,
                     iters=15000, burn=10000, update=500, thresh.all=0,
                     rho.upper=15, nu.upper=10,
                     skew=TRUE, min.s=c(0, 0), max.s=c(10, 10), nknots=5,
                     temporalw=FALSE, temporaltau=FALSE, temporalz=FALSE)
}

# come up with quantile and brier scores
probs <- seq(0.90, 0.99, by=0.01)
quant.scores.gau <- brier.scores.gau <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.1st <- brier.scores.1st <- matrix(NA, nrow=5, ncol=length(probs))
quant.scores.3st <- brier.scores.3st <- matrix(NA, nrow=5, ncol=length(probs))
for (i in 1:5) {
  threshs <- quantile(data[[i]]$y, probs=probs)
  quant.scores.gau[i, ] <- QuantScore(preds=fit.1[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.1st[i, ] <- QuantScore(preds=fit.2[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  quant.scores.3st[i, ] <- QuantScore(preds=fit.3[[i]]$yp, probs=probs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.gau[i, ] <- BrierScore(preds=fit.1[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.1st[i, ] <- BrierScore(preds=fit.2[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
  brier.scores.3st[i, ] <- BrierScore(preds=fit.3[[i]]$yp, thresholds=threshs,
                                      validate=data[[i]]$y[101:144, ])
}

# one bs and qs for each method
quant.score <- brier.score <- matrix(NA, nrow=3, ncol=length(probs))
quant.score[1, ] <- apply(quant.scores.gau, 2, mean)
quant.score[2, ] <- apply(quant.scores.1st, 2, mean)
quant.score[3, ] <- apply(quant.scores.3st, 2, mean)
brier.score[1, ] <- apply(brier.scores.gau, 2, mean)
brier.score[2, ] <- apply(brier.scores.1st, 2, mean)
brier.score[3, ] <- apply(brier.scores.3st, 2, mean)


# Troubleshooting
# Test 1 - Debugging covariance parameters
set.seed(20)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=0, phi.z=0, phi.w=0, phi.tau=0)

source('./mcmc1.R', chdir=T)
source('./auxfunctions.R')
nreps <- 6000
mu <- matrix(10, ns, nt)
res <- data$y - mu
taug <- matrix(NA, ns, nt)
for (t in 1:nt) {
  taug[, t] <- data$tau[, t]
}

rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)

rho <- 5
nu  <- 0.5
gamma <- 0.9
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma.t * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1

par(mfrow=c(1, 2))
set.seed(50)
for (i in 1:nreps) {
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d,
                              gamma=gamma.t, res=res, taug=taug, prec=prec,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh
  }

  if (i %% 200 == 0) {
    plot(rho.keep[1:i], type="l")
    plot(nu.keep[1:i], type="l")
  }
}

gamma <- 0.5
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho.t),
                     error.var=0, smoothness=nu.t, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1
gamma.keep <- rep(NA, nreps)
for (i in 1:nreps) {
  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho.t, nu=nu.t, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  if (i %% 200 == 0) {
    plot(gamma.keep[1:i], type="l")
  }
}

gamma <- 0.5
rho   <- 5
nu    <- 1
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1
rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)
gamma.keep <- rep(NA, nreps)

par(mfrow=c(1, 3))
for (i in 1:nreps) {
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d,
                              gamma=gamma, res=res, taug=taug, prec=prec,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho, nu=nu, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu
  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh

    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  if (i %% 200 == 0) {
    if (i > 5000) {
      start <- i - 5000
    } else {
      start <- 1
    }
    plot(rho.keep[start:i], type="l")
    plot(nu.keep[start:i], type="l")
    plot(gamma.keep[start:i], type="l")
    print(paste("iter", i))
  }
}

source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=1, lambda=0, phi.z=0, phi.w=0, phi.tau=0)
gamma <- 0.5
rho   <- 5
nu    <- 1
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1
nreps <- 10000
rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)
gamma.keep <- rep(NA, nreps)

par(mfrow=c(1, 3))
for (i in 1:nreps) {
  rhonugamma.update <- updateRhoNuGamma(rho=rho, logrho.m=0, logrho.s=10,
                                  fixnu=FALSE, nu=nu, lognu.m=-1.2, lognu.s=1,
                                  d=d, gamma=gamma, res=res, taug=taug,
                                  prec=prec, logdet.prec=logdet.prec,
                                  cur.rss=cur.rss, att.rho=att.rho,
                                  acc.rho=acc.rho, mh.rho=mh.rho,
                                  att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu,
                                  att.gamma=att.gamma, acc.gamma=acc.gamma,
                                  mh.gamma=mh.gamma)
  rho         <- rhonugamma.update$rho
  nu          <- rhonugamma.update$nu
  gamma       <- rhonugamma.update$gamma
  prec        <- rhonugamma.update$prec
  logdet.prec <- rhonugamma.update$logdet.prec
  cur.rss     <- rhonugamma.update$cur.rss
  att.rho     <- rhonugamma.update$att.rho
  acc.rho     <- rhonugamma.update$acc.rho
  att.nu      <- rhonugamma.update$att.nu
  acc.nu      <- rhonugamma.update$acc.nu
  att.gamma   <- rhonugamma.update$att.gamma
  acc.gamma   <- rhonugamma.update$acc.gamma

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu
  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh

    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  if (i %% 200 == 0) {
    if (i > 5000) {
      start <- i - 5000
    } else {
      start <- 1
    }
    plot(rho.keep[start:i], type="l")
    plot(nu.keep[start:i], type="l")
    plot(gamma.keep[start:i], type="l")
    print(paste("iter", i))
  }
}

################################################################################
# Troubleshooting
# Test 4 - Debugging model parameters
################################################################################

# Lambda.1, Lambda.2, and z
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

x.beta   <- matrix(10, ns, nt)
mu       <- x.beta + lambda.1.t * zg

# test lambda.1 updates correctly
lambda.1.keep <- rep(NA, nreps)
for (i in 1:nreps) {
  lambda.1.keep[i] <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                                    taug=taug)
}
plot(lambda.1.keep, type="l")

# test lambda.2 updates correctly
lambda.2.keep <- rep(NA, nreps)
for (i in 1:nreps) {
  lambda.2.keep[i] <- updateLambda2(lambda.a=1, lambda.b=1,
                                    z=data$z, tau=data$tau)
}
plot(lambda.2.keep, type="l", main=lambda.2.t)

# test z updates correctly
z.keep <- array(NA, dim=c(nreps, nknots, nt))
par(mfrow=c(nknots, 5))
for (i in 1:nreps) {
  mu <- x.beta + lambda.1.t * zg
  z.update <- updateZ(y=data$y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=data$tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1.t, lambda.2=lambda.2.t)
  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  if (i %% 500 == 0) {
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (j in 1:nknots) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
             xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    print(paste("iter", i))
  }
}

cover <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover <- cover + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}

# z and lambda updates together
lambda.1.keep <- lambda.2.keep <- rep(NA, nreps)
z.keep <- array(NA, dim=c(nreps, nknots, nt))
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
lambda.1 <- 1
lambda.2 <- 0.5
par(mfrow=c(nknots, 5))
for (i in 1:nreps) {
  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=data$tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg

  z.update <- updateZ(y=data$y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=data$tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)
  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  if (i %% 500 == 0) {
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (j in 1:nknots) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
             xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    print(paste("iter", i))
  }
}

plot(lambda.1.keep, type="l")
plot(lambda.2.keep[1000:6000], type="l", main=lambda.2.t)

cover <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover <- cover + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}

# Beta
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

x.beta   <- matrix(10, ns, nt)

beta.keep <- matrix(NA, nrow=nreps, ncol=3)
par(mfrow=c(1, 3))
for (i in 1:nreps) {
  beta.keep[i, ] <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                               lambda.1=lambda.1.t, taug=taug, prec=prec.t)

  if (i %% 500 == 0) {
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    plot(beta.keep[1:i, 1], type="l")
    plot(beta.keep[1:i, 2], type="l")
    plot(beta.keep[1:i, 3], type="l")
    print(paste("iter", i))
  }
}

# beta, z, and lambda terms
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

# storage
x.beta   <- matrix(10, ns, nt)
beta.keep <- matrix(NA, nrow=nreps, ncol=3)
lambda.1.keep <- lambda.2.keep <- rep(NA, nreps)
z.keep <- array(NA, dim=c(nreps, nknots, nt))

for (i in 1:nreps) {
  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=data$tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg
  z.update <- updateZ(y=data$y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=data$tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)
  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  if (i > 4000) {
    start <- i - 4000
  } else {
    start <- 1
  }

  if ((i + 500) %% 1000 == 0) {
    par(mfrow=c(1, 5))
    plot(beta.keep[start:i, 1], type="l")
    plot(beta.keep[start:i, 2], type="l")
    plot(beta.keep[start:i, 3], type="l")
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
  }
  if (i %% 1000 == 0) {
    par(mfrow=c(nknots, 5))
    for (j in 1:nknots) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
             xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    print(paste("iter", i))
  }
}

par(mfrow=c(1, 3))
plot(beta.keep[, 1], type="l")
bounds <- quantile(beta.keep[, 1], probs=c(0.025, 0.975))
abline(h=bounds)
plot(beta.keep[, 2], type="l")
bounds <- quantile(beta.keep[, 2], probs=c(0.025, 0.975))
abline(h=bounds)
plot(beta.keep[, 3], type="l")
bounds <- quantile(beta.keep[, 3], probs=c(0.025, 0.975))
abline(h=bounds)
# dataset 1: captures beta0
# dataset 2: captures beta0 and beta2, but not beta1

cover <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover <- cover + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}
# dataset 1: coverage around 90%
# dataset 2: coverage around 93%

# Covariance parameters
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

x.beta   <- matrix(10, ns, nt)
mu <- x.beta + lambda.1.t * zg
res <- data$y - mu

rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)

rho <- 5
nu  <- 0.5
gamma <- 0.9
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma.t * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1

par(mfrow=c(1, 2))
set.seed(50)
for (i in 1:nreps) {
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d, rho.upper=15, nu.upper=10,
                              gamma=gamma.t, res=res, taug=taug, prec=prec,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh
  }

  if (i %% 200 == 0) {
    plot(rho.keep[1:i], type="l")
    plot(nu.keep[1:i], type="l")
  }
}

gamma <- 0.5
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho.t),
                     error.var=0, smoothness=nu.t, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1
gamma.keep <- rep(NA, nreps)
for (i in 1:nreps) {
  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho.t, nu=nu.t, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  if (i %% 200 == 0) {
    plot(gamma.keep[1:i], type="l")
  }
}

# Covariance parameters
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

x.beta   <- matrix(10, ns, nt)
mu <- x.beta + lambda.1.t * zg
res <- data$y - mu

gamma <- 0.5
rho   <- 5
nu    <- 1
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1
rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)
gamma.keep <- rep(NA, nreps)

par(mfrow=c(1, 3))
for (i in 1:nreps) {
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d,
                              rho.upper=15, nu.upper=10,
                              gamma=gamma, res=res, taug=taug, prec=prec,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho, nu=nu, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu
  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh

    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  if (i %% 200 == 0) {
    if (i > 5000) {
      start <- i - 5000
    } else {
      start <- 1
    }
    plot(rho.keep[start:i], type="l")
    plot(nu.keep[start:i], type="l")
    plot(gamma.keep[start:i], type="l")
    print(paste("iter", i))
  }
}

# test covariance parameters, beta, and z parameters at one time
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

x.beta   <- matrix(10, ns, nt)
lambda.1 <- 1
mu <- x.beta + lambda.1.t * zg
res <- data$y - mu

gamma <- 0.5
rho   <- 5
nu    <- 1
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1
rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)
gamma.keep <- rep(NA, nreps)
beta.keep  <- matrix(NA, nreps, 3)
lambda.1.keep <- lambda.2.keep <- rep(NA, nreps)
z.keep     <- array(NA, dim=c(nreps, nknots, nt))

for (i in 1:nreps) {
  mu <- x.beta + lambda.1 * zg
  res <- data$y - mu
  cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d,
                              rho.upper=15, nu.upper=10,
                              gamma=gamma, res=res, taug=taug, prec=prec,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  mu <- x.beta + lambda.1 * zg
  res <- data$y - mu
  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho, nu=nu, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu
  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh

    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=data$tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg
  z.update <- updateZ(y=data$y, x.beta=x.beta, zg=zg, prec=prec,
                      tau=data$tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)
  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  if (i %% 200 == 0) {
    if (i > 5000) {
      start <- i - 5000
    } else {
      start <- 1
    }
    par(mfrow=c(2, 3))
    plot(rho.keep[start:i], type="l")
    plot(nu.keep[start:i], type="l")
    plot(gamma.keep[start:i], type="l")
    plot(beta.keep[start:i, 1], type="l")
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    print(paste("iter", i))
  }
}

cover <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover <- cover + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}

# Variance terms
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  zg[, t] <- data$z[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
tau.keep <- array(NA, dim=c(nreps, nknots, nt))

# mh stuff
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)

# initialize testing variables
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
res <- data$y - x.beta - lambda.1.t * zg

for (i in 1:nreps) {
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec.t, z=data$z,
                          lambda.2=lambda.2.t, tau.alpha=tau.alpha.t,
                          tau.beta=tau.beta.t, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau.ns <- tau.update$acc
  att.tau.ns <- tau.update$att
  acc.tau <- tau.update$acc.tau
  att.tau <- tau.update$att.tau



  tau.keep[i, , ] <- tau
  if (i < 4000) {
    start <- 1
  } else {
    start <- i - 4000
  }
  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    for (j in 1:nknots) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(tau.keep[start:i, j, k*10], type="l",
             main=round(data$tau[j, k*10], 3), xlab=nparts)
        bounds <- quantile(tau.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    print(paste("iter", i))
  }
}

cover <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover <- cover + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}

# Variance terms with covariance params
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  zg[, t] <- data$z[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
tau.keep <- array(NA, dim=c(nreps, nknots, nt))
tau.alpha.keep <- tau.beta.keep <- rep(NA, nreps)
rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)
gamma.keep <- rep(NA, nreps)

# mh stuff
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1

# initialize testing variables
tau.alpha <- tau.alpha.t
tau.beta  <- tau.beta.t
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
mu <- x.beta + lambda.1.t * zg
res <- data$y - mu
gamma <- 0.5
rho   <- 5
nu    <- 1
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))

for (i in 1:nreps) {
  tau.alpha <- updateTauAlpha(tau=tau, tau.beta=tau.beta)
  tau.alpha.keep[i] <- tau.alpha

  tau.beta  <- updateTauBeta(tau=tau, tau.alpha=tau.alpha,
                             tau.beta.a=1, tau.beta.b=1)
  tau.beta.keep[i] <- tau.beta

  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec, z=data$z,
                          lambda.2=lambda.2.t, tau.alpha=tau.alpha,
                          tau.beta=tau.beta, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau.ns <- tau.update$acc
  att.tau.ns <- tau.update$att
  acc.tau <- tau.update$acc.tau
  att.tau <- tau.update$att.tau

  tau.keep[i, , ] <- tau

  mu <- x.beta + lambda.1.t * zg
  res <- data$y - mu
  cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d,
                              rho.upper=15, nu.upper=10,
                              gamma=gamma, res=res, taug=taug, prec=prec,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho, nu=nu, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu
  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh

    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  if (i < 4000) {
    start <- 1
  } else {
    start <- i - 4000
  }
  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    plot(tau.alpha.keep[start:i], type="l")
    plot(tau.beta.keep[start:i], type="l")
    plot(rho.keep[start:i], type="l")
    plot(nu.keep[start:i], type="l")
    plot(gamma.keep[start:i], type="l")
    for (j in 1:2) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(tau.keep[start:i, j, k*10], type="l",
             main=round(data$tau[j, k*10], 3), xlab=nparts)
        bounds <- quantile(tau.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    print(paste("iter", i))
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}

# tau.alpha, tau.beta with lambdas
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 2
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- zg <- taug <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  zg[, t] <- data$z[g[, t], t]
  taug[, t] <- data$tau[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
tau.alpha.keep <- tau.beta.keep <- rep(NA, nreps)

# mh stuff
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)

# initialize testing variables
tau.alpha <- 1
tau.beta <- 1
lambda.1 <- 0
lambda.2 <- 1
res <- data$y - x.beta - lambda.1 * zg

for (i in 1:nreps) {
  tau.alpha <- updateTauAlpha(tau=data$tau, tau.beta=tau.beta)
  tau.alpha.keep[i] <- tau.alpha

  tau.beta  <- updateTauBeta(tau=data$tau, tau.alpha=tau.alpha,
                             tau.beta.a=1, tau.beta.b=1)
  tau.beta.keep[i] <- tau.beta

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=data$z, tau=data$tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg

  if (i < 4000) {
    start <- 1
  } else {
    start <- i - 4000
  }
  if (i %% 500 == 0) {
    par(mfrow=c(2, 2))
    plot(tau.alpha.keep[start:i], type="l")
    plot(tau.beta.keep[start:i], type="l")
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    print(paste("iter", i))
  }
}

# Variance terms with lambdas
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 1
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, matrix(data$knots[, , t], 1, 2))
  zg[, t] <- data$z[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
tau.keep <- array(NA, dim=c(nreps, nknots, nt))
tau.alpha.keep <- tau.beta.keep <- rep(NA, nreps)

# mh stuff
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)

# initialize testing variables
tau.alpha <- 1
tau.beta <- 1
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
lambda.1 <- 0
lambda.2 <- 1
res <- data$y - x.beta - lambda.1 * zg

for (i in 1:nreps) {
  tau.alpha <- updateTauAlpha(tau=tau, tau.beta=tau.beta)
  tau.alpha.keep[i] <- tau.alpha

  tau.beta  <- updateTauBeta(tau=tau, tau.alpha=tau.alpha,
                             tau.beta.a=1, tau.beta.b=1)
  tau.beta.keep[i] <- tau.beta

  res <- data$y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec.t, z=data$z,
                          lambda.2=lambda.2, tau.alpha=tau.alpha,
                          tau.beta=tau.beta, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att
  tau.keep[i, , ] <- tau

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=0.1, lambda.b=0.1, z=data$z, tau=tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg

  if (i < 4000) {
    start <- 1
  } else {
    start <- i - 4000
  }
  if (i %% 500 == 0) {
    par(mfrow=c(2, 5))
    for (j in 1) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(tau.keep[start:i, j, k*10], type="l",
             main=round(data$tau[j, k*10], 3), xlab=nparts)
        bounds <- quantile(tau.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    plot(tau.alpha.keep[start:i], type="l")
    plot(tau.beta.keep[start:i], type="l")
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    print(paste("iter", i))
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}

# Variance terms with z
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
tau.keep <- z.keep <- array(NA, dim=c(nreps, nknots, nt))
tau.alpha.keep <- tau.beta.keep <- rep(NA, nreps)

# mh stuff
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)

# initialize testing variables
tau.alpha <- 1
tau.beta <- 1
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
lambda.1 <- 0
lambda.2 <- 1
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
mu <- x.beta + lambda.1 * zg
res <- data$y - mu


for (i in 1:nreps) {
  tau.alpha <- updateTauAlpha(tau=tau, tau.beta=tau.beta)
  tau.alpha.keep[i] <- tau.alpha

  tau.beta  <- updateTauBeta(tau=tau, tau.alpha=tau.alpha,
                             tau.beta.a=1, tau.beta.b=1)
  tau.beta.keep[i] <- tau.beta

  res <- data$y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec.t, z=z,
                          lambda.2=lambda.2, tau.alpha=tau.alpha,
                          tau.beta=tau.beta, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att
  tau.keep[i, , ] <- tau

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg
  z.update <- updateZ(y=data$y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)
  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  if (i < 1000) {
    start <- 1
  } else {
    start <- i - 1000
  }
  if (i %% 500 == 0) {
    par(mfrow=c(3, 5))
    for (j in 1) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(tau.keep[start:i, j, k*10], type="l",
             main=round(data$tau[j, k*10], 3), xlab=nparts)
        bounds <- quantile(tau.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    for (j in 1) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l",
             main=round(data$z[j, k*10], 3), xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    plot(tau.alpha.keep[start:i], type="l")
    plot(tau.beta.keep[start:i], type="l")
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    print(paste("iter", i))
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}

cover.z <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.z <- cover.z + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}

# Variance terms with covariance terms and z
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 10
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
tau.keep <- z.keep <- array(NA, dim=c(nreps, nknots, nt))
tau.alpha.keep <- tau.beta.keep <- rep(NA, nreps)
rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)
gamma.keep <- rep(NA, nreps)

# mh stuff
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1

# initialize testing variables
tau.alpha <- 1
tau.beta <- 1
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
lambda.1 <- 0
lambda.2 <- 1
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
mu <- x.beta + lambda.1 * zg
res <- data$y - mu

gamma <- 0.5
rho   <- 5
nu    <- 1
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))

for (i in 1:nreps) {
  tau.alpha <- updateTauAlpha(tau=tau, tau.beta=tau.beta)
  tau.alpha.keep[i] <- tau.alpha

  tau.beta  <- updateTauBeta(tau=tau, tau.alpha=tau.alpha,
                             tau.beta.a=1, tau.beta.b=1)
  tau.beta.keep[i] <- tau.beta

  mu <- x.beta + lambda.1.t * zg
  res <- data$y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec, z=z,
                          lambda.2=lambda.2, tau.alpha=tau.alpha,
                          tau.beta=tau.beta, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att
  tau.keep[i, , ] <- tau

  mu <- x.beta + lambda.1.t * zg
  res <- data$y - mu
  cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d,
                              gamma=gamma, res=res, taug=taug, prec=prec,
                              rho.upper=15, nu.upper=10,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho, nu=nu, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu
  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh

    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg
  z.update <- updateZ(y=data$y, x.beta=x.beta, zg=zg, prec=prec,
                      tau=tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)
  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  if (i < 4000) {
    start <- 1
  } else {
    start <- i - 4000
  }
  if (i %% 500 == 0) {
    par(mfrow=c(3, 5))
    for (j in 1) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(tau.keep[start:i, j, k*10], type="l",
             main=round(data$tau[j, k*10], 3), xlab=nparts)
        bounds <- quantile(tau.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    for (j in 1) {
      for (k in 1:3) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l",
             main=round(data$z[j, k*10], 3), xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    plot(tau.alpha.keep[start:i], type="l")
    plot(tau.beta.keep[start:i], type="l")
    plot(rho.keep[start:i], type="l")
    plot(nu.keep[start:i], type="l")
    plot(gamma.keep[start:i], type="l")
    print(paste("iter", i))
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}

cover.z <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.z <- cover.z + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}

# Variance terms with covariance terms, z, and beta
set.seed(20)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- -5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 10
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
tau.keep <- z.keep <- array(NA, dim=c(nreps, nknots, nt))
tau.alpha.keep <- tau.beta.keep <- rep(NA, nreps)
rho.keep   <- rep(NA, nreps)
nu.keep    <- rep(NA, nreps)
gamma.keep <- rep(NA, nreps)
beta.keep  <- matrix(NA, nreps, 3)

# mh stuff
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)
att.rho   <- acc.rho   <- mh.rho   <- 0.1
att.nu    <- acc.nu    <- mh.nu    <- 0.1
att.gamma <- acc.gamma <- mh.gamma <- 1

# initialize testing variables
tau.alpha <- 1
tau.beta <- 1
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
lambda.1 <- 0
lambda.2 <- 1
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
mu <- x.beta + lambda.1 * zg
res <- data$y - mu
gamma <- 0.5
rho   <- 5
nu    <- 1
cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                     error.var=0, smoothness=nu, finescale.var=0)
C   <- gamma * cor
diag(C) <- 1
CC <- chol.inv(C, inv=TRUE, logdet=TRUE)
prec <- CC$prec
logdet.prec <- CC$logdet.prec
cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))


for (i in 1:nreps) {
  tau.alpha <- updateTauAlpha(tau=tau, tau.beta=tau.beta)
  tau.alpha.keep[i] <- tau.alpha

  tau.beta  <- updateTauBeta(tau=tau, tau.alpha=tau.alpha,
                             tau.beta.a=1, tau.beta.b=1)
  tau.beta.keep[i] <- tau.beta

  res <- data$y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec, z=z,
                          lambda.2=lambda.2, tau.alpha=tau.alpha,
                          tau.beta=tau.beta, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att
  tau.keep[i, , ] <- tau

  mu <- x.beta + lambda.1 * zg
  res <- data$y - mu
  cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))
  rhonu.update <- updateRhoNu(rho=rho, logrho.m=0, logrho.s=10, fixnu=FALSE,
                              nu=nu, lognu.m=-1.2, lognu.s=1, d=d,
                              rho.upper=15, nu.upper=10,
                              gamma=gamma, res=res, taug=taug, prec=prec,
                              cor=cor, logdet.prec=logdet.prec, cur.rss=cur.rss,
                              att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                              att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
  rho         <- rhonu.update$rho
  nu          <- rhonu.update$nu
  prec        <- rhonu.update$prec
  logdet.prec <- rhonu.update$logdet.prec
  cor         <- rhonu.update$cor
  cur.rss     <- rhonu.update$cur.rss
  att.rho     <- rhonu.update$att.rho
  acc.rho     <- rhonu.update$acc.rho
  att.nu      <- rhonu.update$att.nu
  acc.nu      <- rhonu.update$acc.nu

  gamma.update <- updateGamma(gamma=gamma, gamma.m=0, gamma.s=1,
                              d=d, rho=rho, nu=nu, taug=taug, res=res,
                              prec=prec, cor=cor, logdet.prec=logdet.prec,
                              cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                              mh=mh.gamma)
  gamma       <- gamma.update$gamma
  prec        <- gamma.update$prec
  logdet.prec <- gamma.update$logdet.prec
  acc.gamma   <- gamma.update$acc
  att.gamma   <- gamma.update$att
  cur.rss     <- gamma.update$cur.rss

  rho.keep[i]   <- rho
  nu.keep[i]    <- nu
  gamma.keep[i] <- gamma

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
    acc.nu    <- mh.update$acc
    att.nu    <- mh.update$att
    mh.nu     <- mh.update$mh

    mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
    acc.gamma <- mh.update$acc
    att.gamma <- mh.update$att
    mh.gamma  <- mh.update$mh
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg
  z.update <- updateZ(y=data$y, x.beta=x.beta, zg=zg, prec=prec,
                      tau=tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)
  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  if (i < 4000) {
    start <- 1
  } else {
    start <- i - 4000
  }
  if (i %% 500 == 0) {
    par(mfrow=c(3, 5))
    for (j in 1) {
      for (k in 1:4) {
        nparts <- length(which(g[, k * 10] == j))
        plot(tau.keep[start:i, j, k*10], type="l",
             main=round(data$tau[j, k*10], 3), xlab=nparts)
        bounds <- quantile(tau.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    plot(beta.keep[start:i, 1], type="l")
    for (j in 1) {
      for (k in 1:3) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l",
             main=round(data$z[j, k*10], 3), xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    plot(tau.alpha.keep[start:i], type="l")
    plot(tau.beta.keep[start:i], type="l")
    plot(rho.keep[start:i], type="l")
    plot(nu.keep[start:i], type="l")
    plot(gamma.keep[start:i], type="l")
    print(paste("iter", i))
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}

cover.z <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover.z <- cover.z + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}

################################################################################
# Troubleshooting
# Test 9 - Debugging model parameters
################################################################################
# Lambda.1, Lambda.2, z, and phi.z
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 0.01
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.9
phi.w.t <- 0.8
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

x.beta   <- matrix(10, ns, nt)
mu       <- x.beta + lambda.1.t * zg

# test lambda.1 updates correctly
lambda.1.keep <- rep(NA, nreps)
for (i in 1:nreps) {
  lambda.1.keep[i] <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                                    taug=taug)
}
plot(lambda.1.keep, type="l")

# test lambda.2 updates correctly
lambda.2.keep <- rep(NA, nreps)
for (i in 1:nreps) {
  lambda.2.keep[i] <- updateLambda2(lambda.a=1, lambda.b=1,
                                    z=data$z, tau=data$tau)
}
plot(lambda.2.keep, type="l", main=lambda.2.t)

# test z updates correctly
z.keep <- array(NA, dim=c(nreps, nknots, nt))
z <- matrix(1, nrow=nknots, ncol=nt)
acc.z <- att.z <- mh.z <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
zg <- matrix(1, nrow=ns, ncol=nt)
par(mfrow=c(nknots, 5))
for (i in 1:nreps) {
  mu <- x.beta + lambda.1.t * zg
  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1.t,
                        lambda.2=lambda.2.t, x.beta=x.beta, phi=phi.z.t,
                        tau=data$tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)
  z     <- z.update$z
  zg    <- z.update$zg
  att.z <- z.update$att
  acc.z <- z.update$acc
  z.keep[i, , ] <- z

  if (i %% 500 == 0) {
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (j in 1:nknots) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
             xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z <- mh.update$acc
    att.z <- mh.update$att
    mh.z  <- mh.update$mh
  }
}

cover <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover <- cover + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}
# coverage is around 95%

# test phi.z updates correctly
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.9
phi.w.t <- 0.8
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]    <- data$z[g[, t], t]
}

x.beta   <- matrix(10, ns, nt)
mu       <- x.beta + lambda.1.t * zg
for (t in 1:nt) {
  zg[, t]    <- data$z[g[, t], t]
}
phi.z.keep <- rep(NA, nreps)
phi.z <- 0
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
for (i in 1:nreps) {
  z.update <- updateZTS(z=data$z, zg=zg, y=data$y, lambda.1=lambda.1.t,
                        lambda.2=lambda.2.t, x.beta=x.beta, phi=phi.z,
                        tau=data$tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)

  phi.z     <- z.update$phi
  att.phi.z <- z.update$att
  acc.phi.z <- z.update$acc
  phi.z.keep[i] <- phi.z

  if (i %% 500 == 0) {
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    plot(phi.z.keep[start:i], type="l")
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh
  }
}

# test phi.z, lambda.1, lambda.2, and z updates correctly
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.9
phi.w.t <- 0.8
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
z.keep <- array(NA, dim=c(nreps, nknots, nt))
phi.z.keep <- rep(NA, nreps)

# mh stuff
acc.z <- att.z <- mh.z <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1

# initialize testing variables
z  <- matrix(1, nrow=nknots, ncol=nt)
zg <- matrix(1, nrow=ns, ncol=nt)
mu <- x.beta + lambda.1.t * zg
phi.z <- 0
lambda.1 <- 1
lambda.2 <- 1

for (i in 1:nreps) {
  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=data$tau)

  lambda.2.keep[i] <- lambda.2

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1,
                        lambda.2=lambda.2, x.beta=x.beta, phi=phi.z,
                        tau=data$tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)

  z           <- z.update$z
  zg          <- z.update$zg
  acc.z       <- z.update$acc
  att.z       <- z.update$att
  phi.z       <- z.update$phi
  att.phi.z   <- z.update$att.phi
  acc.phi.z   <- z.update$acc.phi

  # storage
  z.keep[i, , ] <- z
  phi.z.keep[i] <- phi.z

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (j in 1:2) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
             xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z <- mh.update$acc
    att.z <- mh.update$att
    mh.z  <- mh.update$mh
  }
}

cover <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:6000, i, j], probs=c(0.025, 0.975))
    cover <- cover + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}

# test phi.w updates correctly
set.seed(10)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 2
phi.w.t <- 0.8
phi.z.t <- 0.7
phi.tau.t <- 0.7
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]  <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]   <- data$z[g[, t], t]
}
x.beta <- matrix(10, ns, nt)

min.s <- c(0, 0)
max.s <- c(10, 10)

# storage
phi.w.keep <- rep(NA, nreps)

# mh stuff
acc.phi.w <- att.phi.w <- mh.phi.w <- 1
att.w <- acc.w <- mh.w <- matrix(0.15, nknots, nt)

# initialize testing variables
phi.w <- 0
knots.star <- array(NA, dim=c(nknots, 2, nt))
knots.star[, 1, ] <- transform$probit(data$knots[, 1, ], lower=min.s[1],
                                      upper=max.s[1])
knots.star[, 2, ] <- transform$probit(data$knots[, 2, ], lower=min.s[2],
                                      upper=max.s[2])


for (i in 1:nreps) {
  avgparts <- rep(0, nt)
  mu <- x.beta + lambda.1.t * zg
  res <- data$y - mu
  phi.update <- updatePhiTS(data=knots.star, phi=phi.w, day.mar=3,
                            att=att.phi.w, acc=acc.phi.w, mh=mh.phi.w)
  phi.w     <- phi.update$phi
  acc.phi.w <- phi.update$acc.phi
  att.phi.w <- phi.update$att.phi

  phi.w.keep[i] <- phi.w

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.w, att=att.phi.w, mh=mh.phi.w)
    acc.phi.w <- mh.update$acc
    att.phi.w <- mh.update$att
    mh.phi.w  <- mh.update$mh
  }
}
plot(phi.w.keep[1000:nreps], type="l", main=phi.w.t)

# test knots update correctly
set.seed(10)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 2
phi.w.t <- 0.8
phi.z.t <- 0.7
phi.tau.t <- 0.7
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 20000
nparts.tau <- matrix(NA, nknots, nt)
zg <- taug <- g.t <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g.t[, t]  <- mem(s, data$knots[, , t])
}
x.beta <- matrix(10, ns, nt)

# storage
knots.keep <- array(NA, dim=c(nreps, nknots, 2, nt))

# mh stuff
acc.phi.w <- att.phi.w <- mh.phi.w <- 1
att.w <- acc.w <- mh.w <- matrix(0.15, nknots, nt)

# initialize testing variables
min.s <- c(0, 0)
max.s <- c(10, 10)
knots.star <- array(rnorm(nknots * 2 * nt), dim=c(nknots, 2, nt))
knots <- array(NA, dim=c(nknots, 2, nt))
knots[, 1, ] <- transform$inv.probit(knots.star[, 1, ], lower=min.s[1],
                                     upper=max.s[1])
knots[, 2, ] <- transform$inv.probit(knots.star[, 2, ], lower=min.s[2],
                                     upper=max.s[2])
for (t in 1:nt) {  # these need to be initialized with starting knots
  g[, t]    <- mem(s, knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]   <- data$z[g[, t], t]
}

col <- c("firebrick4", "dodgerblue4", "darkolivegreen4")
lcol <- c("firebrick1", "dodgerblue1", "darkolivegreen1")

for (i in 1:nreps) {
  mu <- x.beta + lambda.1.t * zg
  res <- data$y - mu
  knots.update <- updateKnotsTS(phi=phi.w.t, knots=knots, g=g, ts=TRUE,
                                tau=data$tau, z=data$z, s=s, min.s=min.s,
                                max.s=max.s, x.beta=x.beta,
                                lambda.1=lambda.1.t, y=data$y, prec=prec.t,
                                att=att.w, acc=acc.w, mh=mh.w,
                                update.prop=0.25, att.phi=att.phi.w,
                                acc.phi=acc.phi.w, mh.phi=mh.phi.w)
  knots.star <- knots.update$knots.star
  knots      <- knots.update$knots
  g          <- knots.update$g
  taug       <- knots.update$taug
  zg         <- knots.update$zg
  acc.w      <- knots.update$acc
  att.w      <- knots.update$att

  knots.keep[i, , , ] <- knots

  if (i %% 500 == 0) {
    if (i < 1000) {
      start <- 1
    } else {
      start <- i - 1000
    }
    par(mfrow=c(2, 5))
    for (k in 1:5) {
      day <- k * 10
      plot(knots.keep[start:i, 1, , day], type="l", col=lcol[1],
           xlim=c(0, 10), ylim=c(0, 10),
           main=paste("mcmc part, day", day))
      for (j in 2:nknots) {
        lines(knots.keep[start:i, j, , day], col=lcol[j])
      }
      points(knots[, , day], col=col, pch=15, cex=2)
      points(s, col=col[g[, day]])
      cat("\t acc.rate =", round(acc.w[, day] / att.w[, day], 3), "\n")
      cat("\t mh =", round(mh.w[, day], 3), "\n")
    }
    for (k in 1:5) {
      day <- k * 10
      plot(data$knots[, , day], col=col, pch=15, cex=2,
           xlim=c(0, 10), ylim=c(0, 10), main=paste("true part, day", day))
      points(s, xlim=c(0, 10), ylim=c(0, 10), col=col[g.t[, day]])
    }

    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.w, att=att.w, mh=mh.w, nattempts=100)
    acc.w <- mh.update$acc
    att.w <- mh.update$att
    mh.w  <- mh.update$mh
  }

}

classify <- array(0, dim=c(nknots, nknots, nt))
for (t in 1:nt) {
  for (i in 1:ns) {
    row <- g.t[i, t]
    col <- g[i, t]
    classify[row, col, t] <- classify[row, col, t] + 1
  }
}

row1 <- col1 <- c(1, 2, 3)
row2 <- col2 <- c(1, 3, 2)
row3 <- col3 <- c(2, 1, 3)
row4 <- col4 <- c(2, 3, 1)
row5 <- col5 <- c(3, 1, 2)
row6 <- col6 <- c(3, 2, 1)
sums <- matrix(0, nrow=36, ncol=nt)
sum <- rep(0, 36)
for (t in 1:nt) {
  sum[1]  <- sum(diag(classify[row1, col1, t]))
  sum[2]  <- sum(diag(classify[row1, col2, t]))
  sum[3]  <- sum(diag(classify[row1, col3, t]))
  sum[4]  <- sum(diag(classify[row1, col4, t]))
  sum[5]  <- sum(diag(classify[row1, col5, t]))
  sum[6]  <- sum(diag(classify[row1, col6, t]))
  sum[7]  <- sum(diag(classify[row2, col1, t]))
  sum[8]  <- sum(diag(classify[row2, col2, t]))
  sum[9]  <- sum(diag(classify[row2, col3, t]))
  sum[10] <- sum(diag(classify[row2, col4, t]))
  sum[11] <- sum(diag(classify[row2, col5, t]))
  sum[12] <- sum(diag(classify[row2, col6, t]))
  sum[13] <- sum(diag(classify[row3, col1, t]))
  sum[14] <- sum(diag(classify[row3, col2, t]))
  sum[15] <- sum(diag(classify[row3, col3, t]))
  sum[16] <- sum(diag(classify[row3, col4, t]))
  sum[17] <- sum(diag(classify[row3, col5, t]))
  sum[18] <- sum(diag(classify[row3, col6, t]))
  sum[19] <- sum(diag(classify[row4, col1, t]))
  sum[20] <- sum(diag(classify[row4, col2, t]))
  sum[21] <- sum(diag(classify[row4, col3, t]))
  sum[22] <- sum(diag(classify[row4, col4, t]))
  sum[23] <- sum(diag(classify[row4, col5, t]))
  sum[24] <- sum(diag(classify[row4, col6, t]))
  sum[25] <- sum(diag(classify[row5, col1, t]))
  sum[26] <- sum(diag(classify[row5, col2, t]))
  sum[27] <- sum(diag(classify[row5, col3, t]))
  sum[28] <- sum(diag(classify[row5, col4, t]))
  sum[29] <- sum(diag(classify[row5, col5, t]))
  sum[30] <- sum(diag(classify[row5, col6, t]))
  sum[31] <- sum(diag(classify[row6, col1, t]))
  sum[32] <- sum(diag(classify[row6, col2, t]))
  sum[33] <- sum(diag(classify[row6, col3, t]))
  sum[34] <- sum(diag(classify[row6, col4, t]))
  sum[35] <- sum(diag(classify[row6, col5, t]))
  sum[36] <- sum(diag(classify[row6, col6, t]))
  sums[, t] <- sum
}

misclassify <- rep(0, nt)
for (t in 1:nt) {
  misclassify[t] <- 144 - max(sums[, t])
}
sum(misclassify) / (144*50)
# about 12-15% are misclassified over all sites and days when updating each
# knot individually
# about 20% are misclassified over all sites and days when updating all knots on
# a day


# test phi.w and knots updates correctly
set.seed(10)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.w.t <- 0.8
phi.z.t <- 0.7
phi.tau.t <- 0.7
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 7000
nparts.tau <- matrix(NA, nknots, nt)
g.t <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g.t[, t]  <- mem(s, data$knots[, , t])
}
x.beta <- matrix(10, ns, nt)

# storage
phi.w.keep <- rep(NA, nreps)
knots.keep <- array(NA, dim=c(nreps, nknots, 2, nt))

# mh stuff
acc.phi.w <- att.phi.w <- mh.phi.w <- 1
att.w <- acc.w <- mh.w <- matrix(0.15, nknots, nt)

# initialize testing variables
min.s <- c(0, 0)
max.s <- c(10, 10)
fudge <- 1
knots <- data$knots  # this will work alright if we start at the true knots
knots.star <- array(NA, dim=c(nknots, 2, nt))
knots.star[, 1, ] <- transform$probit(x=knots[, 1, ], lower=min.s[1],
                                      upper=max.s[1])
knots.star[, 2, ] <- transform$probit(x=knots[, 2, ], lower=min.s[2],
                                      upper=max.s[2])

knots.star <- array(rnorm(nknots * 2 * nt, knots.star, fudge),
                    dim=c(nknots, 2, nt))
knots[, 1, ] <- transform$inv.probit(x=knots.star[, 1, ], lower=min.s[1],
                                     upper=max.s[1])
knots[, 2, ] <- transform$inv.probit(x=knots.star[, 2, ], lower=min.s[2],
                                     upper=max.s[2])
g <- taug <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {  # these need to be initialized with starting knots
  g[, t]    <- mem(s, knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]   <- data$z[g[, t], t]
}
phi.w <- 0

col <- c("firebrick4", "dodgerblue4", "darkolivegreen4")
lcol <- c("firebrick1", "dodgerblue1", "darkolivegreen1")
for (i in 1:nreps) {
  avgparts <- rep(0, nt)
  knots.update <- updateKnotsTS(phi=phi.w, knots=knots, g=g, ts=TRUE,
                                tau=data$tau, z=data$z, s=s,
                                min.s=min.s, max.s=max.s, x.beta=x.beta,
                                lambda.1=lambda.1.t, y=data$y, prec=prec.t,
                                att=att.w, acc=acc.w, mh=mh.w,
                                update.prop=0.25, att.phi=att.phi.w,
                                acc.phi=acc.phi.w, mh.phi=mh.phi.w)
  knots      <- knots.update$knots
  g          <- knots.update$g
  taug       <- knots.update$taug
  zg         <- knots.update$zg
  acc.z      <- knots.update$acc
  att.z      <- knots.update$att
  phi.w      <- knots.update$phi
  acc.phi.w  <- knots.update$acc.phi
  att.phi.w  <- knots.update$att.phi

  knots.keep[i, , , ] <- knots
  phi.w.keep[i] <- phi.w

  if (i %% 500 == 0) {
    if (i < 1000) {
      start <- 1
    } else {
      start <- i - 1000
    }

    par(mfrow=c(2, 5))
    for (k in 1:4) {
      day <- k * 10
      plot(knots.keep[start:i, 1, , day], type="l", col=lcol[1],
           xlim=c(0, 10), ylim=c(0, 10),
           main=paste("mcmc part, day", day))
      for (j in 2:nknots) {
        lines(knots.keep[start:i, j, , day], col=lcol[j])
      }
      points(knots[, , day], col=col, pch=15, cex=2)
      points(s, col=col[g[, day]])
    }
    plot(phi.w.keep[start:i], type="l", main=phi.w.t)
    for (k in 1:4) {
      day <- k * 10
      plot(data$knots[, , day], col=col, pch=15, cex=2,
           xlim=c(0, 10), ylim=c(0, 10), main=paste("true part, day", day))
      points(s, xlim=c(0, 10), ylim=c(0, 10), col=col[g.t[, day]])
    }

    print(paste("iter", i))

  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.w, att=att.phi.w, mh=mh.phi.w)
    acc.phi.w <- mh.update$acc
    att.phi.w <- mh.update$att
    mh.phi.w  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.w, att=att.w, mh=mh.w, nattempts=100)
    acc.w     <- mh.update$acc
    att.w     <- mh.update$att
    mh.w      <- mh.update$mh
  }
}
plot(phi.w.keep[5001:nreps], type="l")


# test phi.z, z, phi.w, and knots (starting around true knots)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 4
phi.w.t <- 0.8
phi.z.t <- 0.7
phi.tau.t <- 0.7
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 7000
nparts.tau <- matrix(NA, nknots, nt)
g.t <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g.t[, t]  <- mem(s, data$knots[, , t])
}
x.beta <- matrix(10, ns, nt)

# storage
phi.w.keep <- rep(NA, nreps)
phi.z.keep <- rep(NA, nreps)
knots.keep <- array(NA, dim=c(nreps, nknots, 2, nt))
z.keep     <- array(NA, dim=c(nreps, nknots, nt))
zg.keep    <- array(NA, dim=c(nreps, ns, nt))

# mh stuff
acc.phi.w <- att.phi.w <- mh.phi.w <- 1
att.w <- acc.w <- mh.w <- matrix(0.15, nknots, nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
att.z <- acc.z <- mh.z <- matrix(1, nknots, nt)

# initialize testing variables
min.s <- c(0, 0)
max.s <- c(10, 10)
fudge <- 1
knots <- data$knots  # this will work alright if we start at the true knots
knots.star <- array(NA, dim=c(nknots, 2, nt))
knots.star[, 1, ] <- transform$probit(x=knots[, 1, ], lower=min.s[1],
                                      upper=max.s[1])
knots.star[, 2, ] <- transform$probit(x=knots[, 2, ], lower=min.s[2],
                                      upper=max.s[2])

knots.star <- array(rnorm(nknots * 2 * nt, knots.star, fudge),
                    dim=c(nknots, 2, nt))
knots[, 1, ] <- transform$inv.probit(x=knots.star[, 1, ], lower=min.s[1],
                                     upper=max.s[1])
knots[, 2, ] <- transform$inv.probit(x=knots.star[, 2, ], lower=min.s[2],
                                     upper=max.s[2])

z <- matrix(1, nknots, nt)
g <- taug <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {  # these need to be initialized with starting knots
  g[, t]    <- mem(s, knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]   <- z[g[, t], t]
}
phi.w <- 0
phi.z <- 0

col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4")
lcol <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1")
par(mfrow=c(1, 2))
for (i in 1:nreps) {

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1.t,
                        lambda.2=lambda.2.t, x.beta=x.beta, phi=phi.z,
                        tau=data$tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)
  z             <- z.update$z
  zg            <- z.update$zg
  acc.z         <- z.update$acc
  att.z         <- z.update$att
  phi.z         <- z.update$phi
  att.phi.z     <- z.update$att.phi
  acc.phi.z     <- z.update$acc.phi

  phi.z.keep[i] <- phi.z
  z.keep[i, , ] <- z
  zg.keep[i, , ] <- zg

  avgparts <- rep(0, nt)
  knots.update <- updateKnotsTS(phi=phi.w, knots=knots, g=g, ts=TRUE,
                                tau=data$tau, z=z, s=s,
                                min.s=min.s, max.s=max.s,
                                x.beta=x.beta, lambda.1=lambda.1.t, y=data$y,
                                prec=prec.t,  att=att.w, acc=acc.w, mh=mh.w,
                                update.prop=0.25, att.phi=att.phi.w,
                                acc.phi=acc.phi.w, mh.phi=mh.phi.w)
  knots.star <- knots.update$knots.star
  knots      <- knots.update$knots
  g          <- knots.update$g
  taug       <- knots.update$taug
  zg         <- knots.update$zg
  acc.w      <- knots.update$acc
  att.w      <- knots.update$att
  phi.w      <- knots.update$phi
  acc.phi.w  <- knots.update$acc.phi
  att.phi.w  <- knots.update$att.phi

  knots.keep[i, , , ] <- knots
  phi.w.keep[i] <- phi.w

  if (i %% 500 == 0) {
    par(mfrow=c(2, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (k in 1:4) {
      day <- k * 10
      plot(knots.keep[start:i, 1, , day], type="l", col=lcol[1],
           xlim=c(0, 10), ylim=c(0, 10),
           main=paste("mcmc part, day", day))
      for (j in 2:nknots) {
        lines(knots.keep[start:i, j, , day], col=lcol[j])
      }
      points(knots[, , day], col=col, pch=15, cex=2)
      points(s, col=col[g[, day]])
    }
    plot(phi.w.keep[start:i], type="l", main=phi.w.t)
    for (k in 1:4) {
      day <- k * 10
      plot(data$knots[, , day], col=col, pch=15, cex=2,
           xlim=c(0, 10), ylim=c(0, 10), main=paste("true part, day", day))
      points(s, xlim=c(0, 10), ylim=c(0, 10), col=col[g.t[, day]])
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z     <- mh.update$acc
    att.z     <- mh.update$att
    mh.z      <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.w, att=att.w, mh=mh.w, nattempts=nknots * nt)
    acc.w     <- mh.update$acc
    att.w     <- mh.update$att
    mh.w      <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.w, att=att.phi.w, mh=mh.phi.w)
    acc.phi.w <- mh.update$acc
    att.phi.w <- mh.update$att
    mh.phi.w  <- mh.update$mh
  }
}
# updates alright when starting knots are data$knots
# updates alright when starting knots are data$knots with fudge=0.25

g.t <- zg.t <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g.t[, t]    <- mem(s, data$knots[, , t])
  zg.t[, t]   <- data$z[g.t[, t], t]
}

par(mfrow=c(nknots, 5))
for (k in 1:nknots) { for (t in 1:5) {
  site <- k * 30
  day <- t * 10
  plot(zg.keep[, site, day], type="l", main=round(zg.t[site, day], 3))
} }

cover <- 0
for (i in 1:ns) {
  for (t in 1:nt) {
    bounds <- quantile(zg.keep[, i, t], probs=c(0.025, 0.975))
    cover <- cover + ((bounds[1] < zg.t[i, t]) & (bounds[2] > zg.t[i, t]))
  }
}
cover / (ns * nt)
# coverage for zs at sites is around 84% when fudge = 0.5
# coverage for zs at sites is around 80% when fudge = 1.0


# test phi.z, z, phi.w, and knots (starting at random locations)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 4
phi.w.t <- 0.8
phi.z.t <- 0.7
phi.tau.t <- 0.7
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 20000
nparts.tau <- matrix(NA, nknots, nt)
g.t <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g.t[, t]  <- mem(s, data$knots[, , t])
}
x.beta <- matrix(10, ns, nt)

# storage
phi.w.keep <- rep(NA, nreps)
phi.z.keep <- rep(NA, nreps)
knots.keep <- array(NA, dim=c(nreps, nknots, 2, nt))
z.keep     <- array(NA, dim=c(nreps, nknots, nt))
zg.keep    <- array(NA, dim=c(nreps, ns, nt))

# mh stuff
acc.phi.w <- att.phi.w <- mh.phi.w <- 1
att.w <- acc.w <- mh.w <- matrix(0.15, nknots, nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
att.z <- acc.z <- mh.z <- matrix(1, nknots, nt)

# initialize testing variables - trying something new
min.s <- c(0, 0)
max.s <- c(10, 10)
# fudge <- 1
# knots <- data$knots  # this will work alright if we start at the true knots
# knots.star <- array(NA, dim=c(nknots, 2, nt))
# knots.star[, 1, ] <- transform$probit(x=knots[, 1, ], lower=min.s[1],
#                                       upper=max.s[1])
# knots.star[, 2, ] <- transform$probit(x=knots[, 2, ], lower=min.s[2],
#                                       upper=max.s[2])
knots.temp <- array(runif(nknots * 2 * nt, 2, 8), dim=c(nknots, 2, nt))
knots <- array(NA, dim=c(nknots, 2, nt))
g.temp <- taug <- zg <- matrix(NA, ns, nt)

z <- matrix(abs(rnorm(nknots * nt, 0, data$tau)), nknots, nt)
order <- matrix(NA, 24, 4)
order[1, ]  <- c(1, 2, 3, 4)
order[2, ]  <- c(1, 2, 4, 3)
order[3, ]  <- c(1, 3, 2, 4)
order[4, ]  <- c(1, 3, 4, 2)
order[5, ]  <- c(1, 4, 2, 3)
order[6, ]  <- c(1, 4, 3, 2)
order[7, ]  <- c(2, 1, 3, 4)
order[8, ]  <- c(2, 1, 4, 3)
order[9, ]  <- c(2, 3, 1, 4)
order[10, ] <- c(2, 3, 4, 1)
order[11, ] <- c(2, 4, 1, 3)
order[12, ] <- c(2, 4, 3, 1)
order[13, ] <- c(3, 1, 2, 4)
order[14, ] <- c(3, 1, 4, 2)
order[15, ] <- c(3, 2, 1, 4)
order[16, ] <- c(3, 2, 4, 1)
order[17, ] <- c(3, 4, 1, 2)
order[18, ] <- c(3, 4, 2, 1)
order[19, ] <- c(4, 1, 2, 3)
order[20, ] <- c(4, 1, 3, 2)
order[21, ] <- c(4, 2, 1, 3)
order[22, ] <- c(4, 2, 3, 1)
order[23, ] <- c(4, 3, 1, 2)
order[24, ] <- c(4, 3, 2, 1)

for (t in 1:nt) {
  rss <- rep(NA, 6)
  for (o in 1:nrow(order)) {
    g.temp <- mem(s, knots.temp[order[o, ], , t])
    taug.temp <- data$tau[g.temp, t]
    zg.temp   <- z[g.temp, t]
    res.temp  <- data$y[, t] - x.beta[, t] - lambda.1.t * zg.temp
    rss[o] <- quad.form(prec.t, sqrt(taug.temp) * res.temp)
  }
  min.o <- which(rss == min(rss))
  knots[, , t] <- knots.temp[order[min.o, ], , t]
  g[, t] <- mem(s, knots[, , t])
  zg[, t] <- z[g[, t], t]
  taug[, t] <- data$tau[g[, t], t]
}

phi.w <- 0
phi.z <- 0

col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4")
lcol <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1")
par(mfrow=c(1, 2))
for (i in 1:nreps) {

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1.t,
                        lambda.2=lambda.2.t, x.beta=x.beta, phi=phi.z,
                        tau=data$tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)
  z             <- z.update$z
  zg            <- z.update$zg
  acc.z         <- z.update$acc
  att.z         <- z.update$att
  phi.z         <- z.update$phi
  att.phi.z     <- z.update$att.phi
  acc.phi.z     <- z.update$acc.phi

  phi.z.keep[i] <- phi.z
  z.keep[i, , ] <- z
  zg.keep[i, , ] <- zg

  avgparts <- rep(0, nt)
  knots.update <- updateKnotsTS(phi=phi.w, knots=knots, g=g, ts=TRUE,
                                tau=data$tau, z=z, s=s,
                                min.s=min.s, max.s=max.s,
                                x.beta=x.beta, lambda.1=lambda.1.t, y=data$y,
                                prec=prec.t,  att=att.w, acc=acc.w, mh=mh.w,
                                update.prop=0.25, att.phi=att.phi.w,
                                acc.phi=acc.phi.w, mh.phi=mh.phi.w)
  knots.star <- knots.update$knots.star
  knots      <- knots.update$knots
  g          <- knots.update$g
  taug       <- knots.update$taug
  zg         <- knots.update$zg
  acc.w      <- knots.update$acc
  att.w      <- knots.update$att
  phi.w      <- knots.update$phi
  acc.phi.w  <- knots.update$acc.phi
  att.phi.w  <- knots.update$att.phi

  knots.keep[i, , , ] <- knots
  phi.w.keep[i] <- phi.w

  if (i %% 500 == 0) {
    par(mfrow=c(2, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (k in 1:4) {
      day <- k * 10
      plot(knots.keep[start:i, 1, , day], type="l", col=lcol[1],
           xlim=c(0, 10), ylim=c(0, 10),
           main=paste("mcmc part, day", day))
      for (j in 2:nknots) {
        lines(knots.keep[start:i, j, , day], col=lcol[j])
      }
      points(knots[, , day], col=col, pch=15, cex=2)
      points(s, col=col[g[, day]])
    }
    plot(phi.w.keep[start:i], type="l", main=phi.w.t)
    for (k in 1:4) {
      day <- k * 10
      plot(data$knots[, , day], col=col, pch=15, cex=2,
           xlim=c(0, 10), ylim=c(0, 10), main=paste("true part, day", day))
      points(s, xlim=c(0, 10), ylim=c(0, 10), col=col[g.t[, day]])
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    print(paste("iter", i))
  }

  if (i < 15000) {
    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z     <- mh.update$acc
    att.z     <- mh.update$att
    mh.z      <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.w, att=att.w, mh=mh.w, nattempts=nknots * nt)
    acc.w     <- mh.update$acc
    att.w     <- mh.update$att
    mh.w      <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.w, att=att.phi.w, mh=mh.phi.w)
    acc.phi.w <- mh.update$acc
    att.phi.w <- mh.update$att
    mh.phi.w  <- mh.update$mh
  }
}

g.t <- zg.t <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g.t[, t]    <- mem(s, data$knots[, , t])
  zg.t[, t]   <- data$z[g.t[, t], t]
}

par(mfrow=c(nknots, 5))
for (k in 1:nknots) { for (t in 1:5) {
  site <- k * 30
  day <- t * 10
  plot(zg.keep[15001:20000, site, day], type="l",
       main=round(zg.t[site, day], 3))
} }

cover <- 0
for (i in 1:ns) {
  for (t in 1:nt) {
    bounds <- quantile(zg.keep[15001:20000, i, t], probs=c(0.025, 0.975))
    cover <- cover + ((bounds[1] < zg.t[i, t]) & (bounds[2] > zg.t[i, t]))
  }
}
cover / (ns * nt)
# updates alright when starting knots are data$knots
# updates alright when starting knots are data$knots with fudge=0.25
# coverage for zs at sites is around 84% when fudge = 0.5
# coverage for zs at sites is around 80% when fudge = 1.0
# starting with knots labeled on partitions that minimize rss
# with 3 knots coverage for zs at sites is around 78%
# with 4 knots, coverage for zs at sites is around 67%




# test phi.z, z, phi.w, knots, phi.tau, and tau
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.w.t <- 0.7
phi.z.t <- 0.8
phi.tau.t <- 0.8
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 10000
nparts.tau <- matrix(NA, nknots, nt)
g.t <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g.t[, t]  <- mem(s, data$knots[, , t])
}
x.beta <- matrix(10, ns, nt)

# storage
phi.w.keep <- rep(NA, nreps)
phi.z.keep <- rep(NA, nreps)
phi.tau.keep <- rep(NA, nreps)
knots.keep <- array(NA, dim=c(nreps, nknots, 2, nt))
z.keep     <- matrix(NA, nknots, nt)
tau.keep   <- matrix(NA, nknots, nt)

# mh stuff
acc.phi.w <- att.phi.w <- mh.phi.w <- 1
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
acc.phi.tau <- att.phi.tau <- mh.phi.tau <- matrix(NA, nknots, nt)
att.w <- acc.w <- mh.w <- matrix(0.15, nknots, nt)
att.z <- acc.z <- mh.z <- matrix(1, nknots, nt)

# initialize testing variables
min.s <- c(0, 0)
max.s <- c(10, 10)
knots <- array(runif(nknots * 2 * nt, 4, 6), dim=c(nknots, 2, nt))
z     <- matrix(1, nknots, nt)
tau   <- matrix(tau.alpha.t / tau.beta.t, nknots, nt)
g <- taug <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {  # these need to be initialized with starting knots
  g[, t]    <- mem(s, knots[, , t])
  taug[, t] <- tau[g[, t], t]
  zg[, t]   <- z[g[, t], t]
}
phi.z <- 0
phi.w <- 0
phi.tau <- 0

par(mfrow=c(1, 2))
for (i in 1:nreps) {
  mu <- x.beta + lambda.1.t * zg
  tau.update <- res <- data$y - mu
  tau.update <- updateTauTS(phi=phi.tau, tau=tau, taug=taug, g=g, res=res,
                            nparts.tau=nparts.tau, prec=prec.t, z=z,
                            lambda.2=lambda.2.t, tau.alpha=tau.alpha.t,
                            tau.beta=tau.beta.t, skew=TRUE,
                            att=att.tau, acc=acc.tau, mh=mh.tau,
                            att.tau=att.tau, acc.tau=acc.tau,
                            att.phi=att.phi.tau, acc.phi=acc.phi.tau,
                            mh.phi=mh.phi.tau)

  tau     <- tau.update$tau
  taug    <- tau.update$taug
  phi.tau <- tau.update$phi
  acc.tau <- tau.update$acc.tau
  att.tau <- tau.update$att.tau
  att.phi.tau <- tau.update$att.phi
  acc.phi.tau <- tau.update$acc.phi

  tau.keep[i, , ] <- tau
  phi.tau.keep[i] <- phi.tau

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1.t,
                        lambda.2=lambda.2.t, x.beta=x.beta, phi=phi.z,
                        tau=tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)
  z             <- z.update$z
  zg            <- z.update$zg
  acc.z         <- z.update$acc
  att.z         <- z.update$att
  phi.z         <- z.update$phi
  att.phi.z     <- z.update$att.phi
  acc.phi.z     <- z.update$acc.phi

  phi.z.keep[i] <- phi.z

  avgparts <- rep(0, nt)
  knots.update <- updateKnotsTS(phi=phi.w, knots=knots, g=g, ts=TRUE,
                                tau=tau, z=z, s=s,
                                min.s=min.s, max.s=max.s,
                                x.beta=x.beta, lambda.1=lambda.1.t, y=data$y,
                                prec=prec.t,  att=att.w, acc=acc.w, mh=mh.w,
                                update.prop=0.25, att.phi=att.phi.w,
                                acc.phi=acc.phi.w, mh.phi=mh.phi.w)
  knots.star <- knots.update$knots.star
  knots      <- knots.update$knots
  g          <- knots.update$g
  taug       <- knots.update$taug
  zg         <- knots.update$zg
  acc.w      <- knots.update$acc
  att.w      <- knots.update$att
  phi.w      <- knots.update$phi
  acc.phi.w  <- knots.update$acc.phi
  att.phi.w  <- knots.update$att.phi

  knots.keep[i, , , ] <- knots
  phi.w.keep[i] <- phi.w

  if (i %% 500 == 0) {
    par(mfrow=c(3, 4))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (k in 1:4) {
      day <- k * 10
      plot(knots.keep[start:i, 1, , day], type="l", col=lcol[1],
           xlim=c(0, 10), ylim=c(0, 10),
           main=paste("mcmc part, day", day))
      for (j in 2:nknots) {
        lines(knots.keep[start:i, j, , day], col=lcol[j])
      }
      points(knots[, , day], col=col, pch=15, cex=2)
      points(s, col=col[g[, day]])
    }
    for (k in 1:4) {
      day <- k * 10
      plot(data$knots[, , day], col=col, pch=15, cex=2,
           xlim=c(0, 10), ylim=c(0, 10), main=paste("true part, day", day))
      points(s, xlim=c(0, 10), ylim=c(0, 10), col=col[g.t[, day]])
    }
    plot(phi.w.keep[start:i], type="l", main=phi.w.t)
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    plot(phi.tau.keep[start:i], type="l", main=phi.tau.t)
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.tau, att=att.tau, mh=mh.tau)
    acc.tau   <- mh.update$acc
    att.tau   <- mh.update$att
    mh.tau    <- mh.update$mh

    mh.update   <- mhupdate(acc=acc.phi.tau, att=att.phi.tau, mh=mh.phi.tau)
    acc.phi.tau <- mh.update$acc
    att.phi.tau <- mh.update$att
    mh.phi.tau  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z     <- mh.update$acc
    att.z     <- mh.update$att
    mh.z      <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.w, att=att.w, mh=mh.w)
    acc.w     <- mh.update$acc
    att.w     <- mh.update$att
    mh.w      <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.w, att=att.phi.w, mh=mh.phi.w)
    acc.phi.w <- mh.update$acc
    att.phi.w <- mh.update$att
    mh.phi.w  <- mh.update$mh
  }
}

# test phi.z, lambda.1, lambda.2, z, and beta updates correctly
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 5
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.9
phi.w.t <- 0.8
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
taug <- g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
z.keep <- array(NA, dim=c(nreps, nknots, nt))
phi.z.keep <- rep(NA, nreps)
beta.keep  <- matrix(NA, nreps, 3)

# mh stuff
acc.z <- att.z <- mh.z <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1

# initialize testing variables
z  <- matrix(1, nrow=nknots, ncol=nt)
zg <- matrix(1, nrow=ns, ncol=nt)
beta <- c(0, 0, 0)
for (t in 1:nt) {
  x.beta[, t] <- x[, t, ] %*% beta
}
lambda.1 <- 1
mu <- x.beta + lambda.1 * zg
phi.z <- 0
lambda.2 <- 1

for (i in 1:nreps) {
  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=data$tau)

  lambda.2.keep[i] <- lambda.2

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1,
                        lambda.2=lambda.2, x.beta=x.beta, phi=phi.z,
                        tau=data$tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)

  z           <- z.update$z
  zg          <- z.update$zg
  acc.z       <- z.update$acc
  att.z       <- z.update$att
  phi.z       <- z.update$phi
  att.phi.z   <- z.update$att.phi
  acc.phi.z   <- z.update$acc.phi

  # storage
  z.keep[i, , ] <- z
  phi.z.keep[i] <- phi.z

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    for (j in 1:2) {
      for (k in 1:5) {
        nparts <- length(which(g[, k * 10] == j))
        plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
             xlab=nparts)
        bounds <- quantile(z.keep[1:i, j, k*10], probs=c(0.025, 0.975))
        abline(h=bounds)
      }
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    plot(beta.keep[start:i, 1], type="l")
    plot(beta.keep[start:i, 2], type="l")
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z <- mh.update$acc
    att.z <- mh.update$att
    mh.z  <- mh.update$mh
  }
}
# test phi.z, lambda.1, lambda.2, z, and beta updates correctly


# test phi.z, lambda.1, lambda.2, z, beta, and tau with no time series
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.9
phi.w.t <- 0.8
phi.tau.t <- 0
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
z.keep <- tau.keep <- array(NA, dim=c(nreps, nknots, nt))
phi.z.keep <- rep(NA, nreps)
beta.keep  <- matrix(NA, nreps, 3)

# mh stuff
acc.z <- att.z <- mh.z <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)

# initialize testing variables
z  <- matrix(1, nrow=nknots, ncol=nt)
zg <- matrix(1, nrow=ns, ncol=nt)
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
beta <- c(0, 0, 0)
for (t in 1:nt) {
  x.beta[, t] <- x[, t, ] %*% beta
}
phi.z <- 0
lambda.2 <- 1

for (i in 1:nreps) {
  mu <- x.beta + lambda.1 * zg
  res <- data$y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec.t, z=z,
                          lambda.2=lambda.2.t, tau.alpha=tau.alpha.t,
                          tau.beta=tau.beta.t, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau.ns <- tau.update$acc
  att.tau.ns <- tau.update$att
  acc.tau <- tau.update$acc.tau
  att.tau <- tau.update$att.tau
  tau.keep[i, , ] <- tau


  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)

  lambda.2.keep[i] <- lambda.2

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1,
                        lambda.2=lambda.2, x.beta=x.beta, phi=phi.z,
                        tau=tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)

  z           <- z.update$z
  zg          <- z.update$zg
  acc.z       <- z.update$acc
  att.z       <- z.update$att
  phi.z       <- z.update$phi
  att.phi.z   <- z.update$att.phi
  acc.phi.z   <- z.update$acc.phi

  # storage
  z.keep[i, , ] <- z
  phi.z.keep[i] <- phi.z

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    j <- 1
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
           xlab=nparts)
      bounds <- quantile(z.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(tau.keep[start:i, j, k*10], type="l",
           main=round(data$tau[j, k*10], 3), xlab=nparts)
      bounds <- quantile(tau.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    plot(beta.keep[start:i, 1], type="l")
    plot(beta.keep[start:i, 2], type="l")
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z <- mh.update$acc
    att.z <- mh.update$att
    mh.z  <- mh.update$mh
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}
# tau coverage is around 69%

cover.z <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.z <- cover.z + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}
# z coverage is around 99%


# test phi.z, lambda.1, lambda.2, z, beta, and no time series on tau
# trying to use the update function though for the tau time series
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.9
phi.w.t <- 0.8
phi.tau.t <- 0
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
z.keep <- tau.keep <- array(NA, dim=c(nreps, nknots, nt))
phi.z.keep <- rep(NA, nreps)
beta.keep  <- matrix(NA, nreps, 3)

# mh stuff
acc.z <- att.z <- mh.z <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.tau <- att.phi.tau <- mh.phi.tau <- 1

# initialize testing variables
z  <- matrix(1, nrow=nknots, ncol=nt)
zg <- matrix(1, nrow=ns, ncol=nt)
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
beta <- c(0, 0, 0)
for (t in 1:nt) {
  x.beta[, t] <- x[, t, ] %*% beta
}
phi.z <- 0
lambda.2 <- 1

for (i in 1:nreps) {
  mu <- x.beta + lambda.1 * zg
  res <- data$y - mu
  tau.update <- updateTauTS(phi=phi.tau.t, tau=tau, taug=taug, g=g, res=res,
                            nparts.tau=nparts.tau, prec=prec.t, z=z,
                            lambda.2=lambda.2.t, tau.alpha=tau.alpha.t,
                            tau.beta=tau.beta.t, skew=TRUE,
                            att=att.tau, acc=acc.tau, mh=mh.tau,
                            att.phi=att.phi.tau, acc.phi=acc.phi.tau,
                            mh.phi=mh.phi.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att
  tau.keep[i, , ] <- tau


  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)

  lambda.2.keep[i] <- lambda.2

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1,
                        lambda.2=lambda.2, x.beta=x.beta, phi=phi.z,
                        tau=tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)

  z           <- z.update$z
  zg          <- z.update$zg
  acc.z       <- z.update$acc
  att.z       <- z.update$att
  phi.z       <- z.update$phi
  att.phi.z   <- z.update$att.phi
  acc.phi.z   <- z.update$acc.phi

  # storage
  z.keep[i, , ] <- z
  phi.z.keep[i] <- phi.z

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    j <- 1
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
           xlab=nparts)
      bounds <- quantile(z.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(tau.keep[start:i, j, k*10], type="l",
           main=round(data$tau[j, k*10], 3), xlab=nparts)
      bounds <- quantile(tau.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    plot(lambda.1.keep[start:i], type="l")
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    plot(beta.keep[start:i, 1], type="l")
    plot(beta.keep[start:i, 2], type="l")
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z <- mh.update$acc
    att.z <- mh.update$att
    mh.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.tau, att=att.phi.tau, mh=mh.phi.tau)
    acc.phi.tau <- mh.update$acc
    att.phi.tau <- mh.update$att
    mh.phi.tau  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.tau, att=att.tau, mh=mh.tau)
    acc.tau <- mh.update$acc
    att.tau <- mh.update$att
    mh.tau  <- mh.update$mh
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}
# coverage is at 69%

cover.z <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.z <- cover.z + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}
# coverage is at 99%


# test phi.z, lambda.1, lambda.2, z, beta, and no time series on tau
# trying to use the update function though for the tau time series
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
set.seed(10)
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.7
phi.w.t <- 0.8
phi.tau.t <- 0
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
z.keep <- tau.keep <- array(NA, dim=c(nreps, nknots, nt))
phi.z.keep <- rep(NA, nreps)
beta.keep  <- matrix(NA, nreps, 3)

# mh stuff
acc.z <- att.z <- mh.z <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.tau <- att.phi.tau <- mh.phi.tau <- 1

# initialize testing variables
z  <- matrix(1, nrow=nknots, ncol=nt)
zg <- matrix(1, nrow=ns, ncol=nt)
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
beta <- c(0, 0, 0)
for (t in 1:nt) {
  x.beta[, t] <- x[, t, ] %*% beta
}
phi.z <- 0
lambda.2 <- 1

for (i in 1:nreps) {
  mu <- x.beta + lambda.1 * zg
  res <- data$y - mu
  tau.update <- updateTauTS(phi=phi.tau.t, tau=tau, taug=taug, g=g, res=res,
                            nparts.tau=nparts.tau, prec=prec.t, z=z,
                            lambda.2=lambda.2, tau.alpha=tau.alpha.t,
                            tau.beta=tau.beta.t, skew=TRUE,
                            att=att.tau, acc=acc.tau, mh=mh.tau,
                            att.phi=att.phi.tau, acc.phi=acc.phi.tau,
                            mh.phi=mh.phi.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att

  # storage
  tau.keep[i, , ] <- tau

  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)

  lambda.2.keep[i] <- lambda.2

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1,
                        lambda.2=lambda.2, x.beta=x.beta, phi=phi.z,
                        tau=tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)

  z           <- z.update$z
  zg          <- z.update$zg
  acc.z       <- z.update$acc
  att.z       <- z.update$att
  phi.z       <- z.update$phi
  att.phi.z   <- z.update$att.phi
  acc.phi.z   <- z.update$acc.phi

  # storage
  z.keep[i, , ] <- z
  phi.z.keep[i] <- phi.z

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    j <- 1
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
           xlab=nparts)
      bounds <- quantile(z.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(tau.keep[start:i, j, k*10], type="l",
           main=round(data$tau[j, k*10], 3), xlab=nparts)
      bounds <- quantile(tau.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    plot(lambda.1.keep[start:i] / sqrt(lambda.2.keep[start:i]), type="l",
         main="lambda")
    plot(beta.keep[start:i, 1], type="l")
    plot(beta.keep[start:i, 2], type="l")
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z <- mh.update$acc
    att.z <- mh.update$att
    mh.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.tau, att=att.tau, mh=mh.tau)
    acc.tau <- mh.update$acc
    att.tau <- mh.update$att
    mh.tau  <- mh.update$mh
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}
# coverage is at 68%

cover.z <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.z <- cover.z + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}
# coverage is at 97%


# test phi.z, lambda.1, lambda.2, z, beta, phi.tau, tau
# trying to use the update function though for the tau time series
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
set.seed(30)
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
phi.z.t <- 0.7
phi.w.t <- 0.8
phi.tau.t <- 0.7
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t] <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
lambda.1.keep <- rep(NA, nreps)
lambda.2.keep <- rep(NA, nreps)
z.keep <- tau.keep <- array(NA, dim=c(nreps, nknots, nt))
phi.z.keep <- phi.tau.keep <- rep(NA, nreps)
beta.keep  <- matrix(NA, nreps, 3)

# mh stuff
acc.z <- att.z <- mh.z <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.z <- att.phi.z <- mh.phi.z <- 1
acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)
acc.phi.tau <- att.phi.tau <- mh.phi.tau <- 1

# initialize testing variables
z  <- matrix(1, nrow=nknots, ncol=nt)
zg <- matrix(1, nrow=ns, ncol=nt)
tau <- matrix(0.375, nrow=nknots, ncol=nt)
taug <- matrix(0.375, nrow=ns, ncol=nt)
beta <- c(0, 0, 0)
for (t in 1:nt) {
  x.beta[, t] <- x[, t, ] %*% beta
}
phi.z <- 0
phi.tau <- 0
lambda.1 <- 1
lambda.2 <- 1

for (i in 1:nreps) {
  mu <- x.beta + lambda.1 * zg
  res <- data$y - mu
  tau.update <- updateTauTS(phi=phi.tau, tau=tau, taug=taug, g=g, res=res,
                            nparts.tau=nparts.tau, prec=prec.t, z=z,
                            lambda.2=lambda.2, tau.alpha=tau.alpha.t,
                            tau.beta=tau.beta.t, skew=TRUE,
                            att=att.tau, acc=acc.tau, mh=mh.tau,
                            att.phi=att.phi.tau, acc.phi=acc.phi.tau,
                            mh.phi=mh.phi.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att
  phi.tau <- tau.update$phi
  acc.phi.tau <- tau.update$acc.phi
  att.phi.tau <- tau.update$att.phi

  # storage
  tau.keep[i, , ] <- tau
  phi.tau.keep[i] <- phi.tau

  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=data$y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=data$y, prec=prec.t,
                            taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)

  lambda.2.keep[i] <- lambda.2

  z.update <- updateZTS(z=z, zg=zg, y=data$y, lambda.1=lambda.1,
                        lambda.2=lambda.2, x.beta=x.beta, phi=phi.z,
                        tau=tau, taug=taug, g=g, prec=prec.t,
                        acc=acc.z, att=att.z, mh=mh.z,
                        acc.phi=acc.phi.z, att.phi=att.phi.z, mh.phi=mh.phi.z)

  z           <- z.update$z
  zg          <- z.update$zg
  acc.z       <- z.update$acc
  att.z       <- z.update$att
  phi.z       <- z.update$phi
  att.phi.z   <- z.update$att.phi
  acc.phi.z   <- z.update$acc.phi

  # storage
  z.keep[i, , ] <- z
  phi.z.keep[i] <- phi.z

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    j <- 1
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(z.keep[start:i, j, k*10], type="l", main=round(data$z[j, k*10], 3),
           xlab=nparts)
      bounds <- quantile(z.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    for (k in 1:5) {
      nparts <- length(which(g[, k * 10] == j))
      plot(tau.keep[start:i, j, k*10], type="l",
           main=round(data$tau[j, k*10], 3), xlab=nparts)
      bounds <- quantile(tau.keep[start:i, j, k*10], probs=c(0.025, 0.975))
      abline(h=bounds)
    }
    plot(phi.z.keep[start:i], type="l", main=phi.z.t)
    plot(phi.tau.keep[start:i], type="l", main=phi.tau.t)
    plot(lambda.1.keep[start:i] / sqrt(lambda.2.keep[start:i]), type="l",
         main="lambda")
    plot(beta.keep[start:i, 1], type="l")
    plot(beta.keep[start:i, 2], type="l")
    print(paste("iter", i))
  }

  if (i < 5000) {
    mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
    acc.phi.z <- mh.update$acc
    att.phi.z <- mh.update$att
    mh.phi.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
    acc.z <- mh.update$acc
    att.z <- mh.update$att
    mh.z  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.phi.tau, att=att.phi.tau, mh=mh.phi.tau)
    acc.phi.tau <- mh.update$acc
    att.phi.tau <- mh.update$att
    mh.phi.tau  <- mh.update$mh

    mh.update <- mhupdate(acc=acc.tau, att=att.tau, mh=mh.tau)
    acc.tau <- mh.update$acc
    att.tau <- mh.update$att
    mh.tau  <- mh.update$mh
  }
}

cover.tau <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.tau <- data$tau[i, j]
    bounds <- quantile(tau.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.tau <- cover.tau + ((this.tau > bounds[1]) & (this.tau < bounds[2])) / (nknots * nt)
  }
}
# seed(10) coverage is at 61%
# seed(20) coverage is at 69%

cover.z <- 0
for (i in 1:nknots) {
  for (j in 1:nt) {
    this.z <- data$z[i, j]
    bounds <- quantile(z.keep[5001:nreps, i, j], probs=c(0.025, 0.975))
    cover.z <- cover.z + ((this.z > bounds[1]) & (this.z < bounds[2])) / (nknots * nt)
  }
}
# seed(10) coverage is at 99%
# seed(20) coverage is at 98%

# checking data imputation
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- taug <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]    <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]   <- data$z[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)
mu <- x.beta + lambda.1.t * zg

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))

# initialize testing variables
thresh.mtx <- matrix(quantile(data$y, probs=0.80), ns, nt)
impute.obs <- data$y <= thresh.mtx
y <- data$y
y[impute.obs] <- thresh.mtx[impute.obs]

pct.below <- apply(impute.obs, 2, mean)

for (i in 1:nreps) {
  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, cor=cor.t, gamma=gamma.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 100 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 4))
    plot(y.keep[start:i, 3, 1], type="l", main="y 1, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[3, 1], 2))
    plot(y.keep[start:i, 5, 2], type="l", main="y 5, 2",
         xlab=round(pct.below[2], 3), ylab=round(data$y[5, 2], 2))
    plot(y.keep[start:i, 5, 8], type="l", main="y 5, 4",
         xlab=round(pct.below[8], 3), ylab=round(data$y[5, 8], 2))
    plot(y.keep[start:i, 5, 45], type="l", main="y 5, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[5, 45], 2))
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 2], data$y[, 2])
    hist.3.range <- range(y.keep[i, , 8], data$y[, 8])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 2], main="imputed day 2", xlim=hist.2.range)
    hist(y.keep[i, , 8], main="imputed day 8", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 2], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 8], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# checking data imputation with beta
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- taug <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]    <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
  zg[, t]   <- data$z[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)
mu <- x.beta + lambda.1.t * zg

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))
beta.keep <- matrix(NA, nreps, 3)

# initialize testing variables
thresh.mtx <- matrix(quantile(data$y, probs=0.80), ns, nt)
y <- ifelse(
  data$y < thresh.mtx,
  thresh.mtx,
  data$y
)
impute.obs <- y <= thresh.mtx
pct.below <- apply(impute.obs, 2, mean)

for (i in 1:nreps) {
  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=y, zg=zg,
                     lambda.1=lambda.1.t, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }
  mu <- x.beta + lambda.1.t * zg

  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, cor=cor.t, gamma=gamma.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 4))
    plot(y.keep[start:i, 3, 1], type="l", main="y 1, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[3, 1], 2))
    plot(y.keep[start:i, 5, 2], type="l", main="y 5, 2",
         xlab=round(pct.below[2], 3), ylab=round(data$y[5, 2], 2))
    plot(y.keep[start:i, 5, 8], type="l", main="y 5, 4",
         xlab=round(pct.below[8], 3), ylab=round(data$y[5, 8], 2))
    plot(y.keep[start:i, 5, 45], type="l", main="y 5, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[5, 45], 2))
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 2], data$y[, 2])
    hist.3.range <- range(y.keep[i, , 8], data$y[, 8])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 2], main="imputed day 2", xlim=hist.2.range)
    hist(y.keep[i, , 8], main="imputed day 8", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 2], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 8], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# checking data imputation with z
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- taug <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]    <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))
z.keep <- array(NA, dim=c(nreps, nknots, nt))

# initialize testing variables
thresh.mtx <- matrix(quantile(data$y, probs=0.80), ns, nt)
y <- ifelse(
  data$y < thresh.mtx,
  thresh.mtx,
  data$y
)
impute.obs <- y <= thresh.mtx
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
mu <- x.beta + lambda.1.t * zg
pct.below <- apply(impute.obs, 2, mean)

for (i in 1:nreps) {
  z.update <- updateZ(y=y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=data$tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1.t, lambda.2=lambda.2.t)

  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  mu <- x.beta + lambda.1.t * zg

  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, cor=cor.t, gamma=gamma.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 4))
    plot(y.keep[start:i, 3, 1], type="l", main="y 1, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[3, 1], 2))
    plot(y.keep[start:i, 5, 2], type="l", main="y 5, 2",
         xlab=round(pct.below[2], 3), ylab=round(data$y[5, 2], 2))
    plot(y.keep[start:i, 5, 8], type="l", main="y 5, 4",
         xlab=round(pct.below[8], 3), ylab=round(data$y[5, 8], 2))
    plot(y.keep[start:i, 5, 45], type="l", main="y 5, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[5, 45], 2))
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 2], data$y[, 2])
    hist.3.range <- range(y.keep[i, , 8], data$y[, 8])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 2], main="imputed day 2", xlim=hist.2.range)
    hist(y.keep[i, , 8], main="imputed day 8", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 2], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 8], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# checking data imputation with z and lambda
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- taug <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]    <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))
z.keep <- array(NA, dim=c(nreps, nknots, nt))
beta.keep <- matrix(NA, nreps, 3)
lambda.1.keep <- lambda.2.keep <- rep(NA, nreps)

# initialize testing variables
y <- data$y
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
mu <- x.beta + lambda.1.t * zg
thresh.mtx <- matrix(quantile(y, probs=0.80), ns, nt)
impute.obs <- y <= thresh.mtx
pct.below <- apply(impute.obs, 2, mean)

for (i in 1:nreps) {
  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=y, prec=prec.t, taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=data$tau)
  lambda.2.keep[i] <- lambda.2

  z.update <- updateZ(y=y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=data$tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)

  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  mu <- x.beta + lambda.1 * zg

  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, cor=cor.t, gamma=gamma.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 4))
    plot(y.keep[start:i, 3, 1], type="l", main="y 1, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[3, 1], 2))
    plot(y.keep[start:i, 5, 2], type="l", main="y 5, 2",
         xlab=round(pct.below[2], 3), ylab=round(data$y[5, 2], 2))
    plot(y.keep[start:i, 5, 8], type="l", main="y 5, 4",
         xlab=round(pct.below[8], 3), ylab=round(data$y[5, 8], 2))
    plot(y.keep[start:i, 5, 45], type="l", main="y 5, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[5, 45], 2))
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 2], data$y[, 2])
    hist.3.range <- range(y.keep[i, , 8], data$y[, 8])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 2], main="imputed day 2", xlim=hist.2.range)
    hist(y.keep[i, , 8], main="imputed day 8", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 2], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 8], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# checking data imputation with z, lambda, and beta
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- taug <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]    <- mem(s, data$knots[, , t])
  taug[, t] <- data$tau[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))
z.keep <- array(NA, dim=c(nreps, nknots, nt))
beta.keep <- matrix(NA, nreps, 3)
lambda.1.keep <- lambda.2.keep <- rep(NA, nreps)

# initialize testing variables
thresh.mtx <- matrix(quantile(data$y, probs=0.80), ns, nt)
y <- ifelse(
  data$y < thresh.mtx,
  thresh.mtx,
  data$y
)
impute.obs <- y <= thresh.mtx
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
mu <- x.beta + lambda.1.t * zg
pct.below <- apply(impute.obs, 2, mean)

for (i in 1:nreps) {
  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=y, prec=prec.t, taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=data$tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg

  z.update <- updateZ(y=y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=data$tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)

  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  mu <- x.beta + lambda.1 * zg

  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, cor=cor.t, gamma=gamma.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 5))
    plot(y.keep[start:i, 3, 1], type="l", main="y 1, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[3, 1], 2))
    plot(y.keep[start:i, 5, 2], type="l", main="y 5, 2",
         xlab=round(pct.below[2], 3), ylab=round(data$y[5, 2], 2))
    plot(y.keep[start:i, 5, 8], type="l", main="y 5, 4",
         xlab=round(pct.below[8], 3), ylab=round(data$y[5, 8], 2))
    plot(y.keep[start:i, 5, 45], type="l", main="y 5, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[5, 45], 2))
    plot(beta.keep[start: i, 1], type="l", main="beta")
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 2], data$y[, 2])
    hist.3.range <- range(y.keep[i, , 8], data$y[, 8])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 2], main="imputed day 2", xlim=hist.2.range)
    hist(y.keep[i, , 8], main="imputed day 8", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 2], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 8], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# checking data imputation with tau
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 1
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
data$z <- matrix(data$z, 1, nt)
g <- zg <- matrix(NA, ns, nt)

for (t in 1:nt) {
  g[, t]    <- 1
  zg[, t]   <- data$z[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)
mu <- x.beta + lambda.1.t * zg

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))
tau.keep <- array(NA, dim=c(nreps, nknots, nt))

# mh adjustments
att.tau <- acc.tau <- mh.tau <- matrix(1, nknots, nt)

# initialize testing variables
thresh.mtx <- matrix(quantile(data$y, probs=0.70), ns, nt)
y <- ifelse(
  data$y < thresh.mtx,
  thresh.mtx,
  data$y
)
impute.obs <- y <= thresh.mtx
pct.below <- apply(impute.obs, 2, mean)
tau <- matrix(0.375, nknots, nt)
taug <- matrix(0.375, ns, nt)

for (i in 1:nreps) {
  mu <- x.beta + lambda.1.t * zg
  res <- y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec.t, z=data$z,
                          lambda.2=lambda.2.t, tau.alpha=tau.alpha.t,
                          tau.beta=tau.beta.t, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att

  tau.keep[i, , ] <- tau

  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, cor=cor.t, gamma=gamma.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 50 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 4))
    plot(y.keep[start:i, 30, 1], type="l", main="y 30, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[30, 1], 2))
    plot(y.keep[start:i, 3, 7], type="l", main="y 3, 7",
         xlab=round(pct.below[7], 3), ylab=round(data$y[3, 7], 2))
    plot(y.keep[start:i, 2, 39], type="l", main="y 2, 39",
         xlab=round(pct.below[39], 3), ylab=round(data$y[2, 39], 2))
    plot(y.keep[start:i, 4, 45], type="l", main="y 4, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[4, 45], 2))
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 7], data$y[, 7])
    hist.3.range <- range(y.keep[i, , 39], data$y[, 39])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 7], main="imputed day 7", xlim=hist.2.range)
    hist(y.keep[i, , 39], main="imputed day 39", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 7], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 39], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# checking data imputation with tau (3 knots)
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- zg <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]    <- mem(s, data$knots[, , t])
  zg[, t]   <- data$z[g[, t], t]
}
x.beta   <- matrix(10, ns, nt)
mu <- x.beta + lambda.1.t * zg

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))
tau.keep <- array(NA, dim=c(nreps, nknots, nt))

# mh adjustments
att.tau <- acc.tau <- mh.tau <- matrix(1, nknots, nt)

# initialize testing variables
thresh.mtx <- matrix(quantile(data$y, probs=0.30), ns, nt)
y <- ifelse(
  data$y < thresh.mtx,
  thresh.mtx,
  data$y
)
impute.obs <- y <= thresh.mtx
pct.below <- apply(impute.obs, 2, mean)
tau <- matrix(0.375, nknots, nt)
taug <- matrix(0.375, ns, nt)

for (i in 1:nreps) {
  mu <- x.beta + lambda.1.t * zg
  res <- y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec.t, z=data$z,
                          lambda.2=lambda.2.t, tau.alpha=tau.alpha.t,
                          tau.beta=tau.beta.t, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att

  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, cor=cor.t, gamma=gamma.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 50 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 4))
    plot(y.keep[start:i, 30, 1], type="l", main="y 30, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[30, 1], 2))
    plot(y.keep[start:i, 3, 7], type="l", main="y 3, 7",
         xlab=round(pct.below[7], 3), ylab=round(data$y[3, 7], 2))
    plot(y.keep[start:i, 2, 39], type="l", main="y 2, 39",
         xlab=round(pct.below[39], 3), ylab=round(data$y[2, 39], 2))
    plot(y.keep[start:i, 4, 45], type="l", main="y 4, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[4, 45], 2))
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 7], data$y[, 7])
    hist.3.range <- range(y.keep[i, , 39], data$y[, 39])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 7], main="imputed day 7", xlim=hist.2.range)
    hist(y.keep[i, , 39], main="imputed day 39", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 7], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 39], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# checking data imputation with z, lambda, beta, and tau
source('./mcmc.R', chdir=T)
source('./auxfunctions.R')
lambda <- 2
lambda.1.t <- sign(lambda)
lambda.2.t <- 1 / lambda^2
nknots <- 3
set.seed(25)
data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t, gamma=gamma.t, nu=nu.t,
                   rho=rho.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                   dist="t", nknots=nknots, lambda=lambda,
                   phi.z=0, phi.w=0, phi.tau=0)

nreps <- 6000
nparts.tau <- matrix(NA, nknots, nt)
g <- matrix(NA, ns, nt)
for (t in 1:nt) {
  g[, t]    <- mem(s, data$knots[, , t])
}
x.beta   <- matrix(10, ns, nt)

# storage
y.keep <- array(NA, dim=c(nreps, ns, nt))
z.keep <- array(NA, dim=c(nreps, nknots, nt))
tau.keep <- array(NA, dim=c(nreps, nknots, nt))
beta.keep <- matrix(NA, nreps, 3)
lambda.1.keep <- lambda.2.keep <- rep(NA, nreps)

# mh adjustments
att.tau <- acc.tau <- mh.tau <- matrix(1, nknots, nt)

# initialize testing variables
thresh.mtx <- matrix(quantile(data$y, probs=0.20), ns, nt)
y <- ifelse(
  data$y < thresh.mtx,
  thresh.mtx,
  data$y
)
impute.obs <- y <= thresh.mtx
z <- matrix(1, nknots, nt)
zg <- matrix(1, ns, nt)
tau <- matrix(0.375, nknots, nt)
taug <- matrix(0.375, ns, nt)
mu <- x.beta + lambda.1.t * zg
pct.below <- apply(impute.obs, 2, mean)

for (i in 1:nreps) {
  mu <- x.beta + lambda.1 * zg
  res <- y - mu
  tau.update <- updateTau(tau=tau, taug=taug, g=g, res=res,
                          nparts.tau=nparts.tau, prec=prec.t, z=z,
                          lambda.2=lambda.2, tau.alpha=tau.alpha.t,
                          tau.beta=tau.beta.t, skew=TRUE,
                          att=att.tau, acc=acc.tau, mh=mh.tau)
  tau  <- tau.update$tau
  taug <- tau.update$taug
  acc.tau <- tau.update$acc
  att.tau <- tau.update$att

  beta <- updateBeta(beta.m=0, beta.s=10, x=x, y=y, zg=zg,
                     lambda.1=lambda.1, taug=taug, prec=prec.t)
  beta.keep[i, ] <- beta

  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=y, prec=prec.t, taug=taug)
  lambda.1.keep[i] <- lambda.1

  lambda.2 <- updateLambda2(lambda.a=1, lambda.b=1, z=z, tau=tau)
  lambda.2.keep[i] <- lambda.2

  mu <- x.beta + lambda.1 * zg

  z.update <- updateZ(y=y, x.beta=x.beta, zg=zg, prec=prec.t,
                      tau=tau, mu=mu, taug=taug, g=g,
                      lambda.1=lambda.1, lambda.2=lambda.2)

  z <- z.update$z
  zg <- z.update$zg
  z.keep[i, , ] <- z

  mu <- x.beta + lambda.1 * zg

  y <- imputeY(y=y, taug=taug, mu=mu, obs=impute.obs, prec=prec.t,
               thresh.mtx=thresh.mtx)
  y.keep[i, , ] <- y

  if (i %% 500 == 0) {
    par(mfrow=c(nknots, 5))
    if (i > 4000) {
      start <- i - 4000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 5))
    plot(y.keep[start:i, 30, 1], type="l", main="y 1, 1",
         xlab=round(pct.below[1], 3), ylab=round(data$y[30, 1], 2))
    plot(y.keep[start:i, 3, 7], type="l", main="y 5, 2",
         xlab=round(pct.below[7], 3), ylab=round(data$y[3, 7], 2))
    plot(y.keep[start:i, 2, 39], type="l", main="y 5, 4",
         xlab=round(pct.below[39], 3), ylab=round(data$y[2, 39], 2))
    plot(y.keep[start:i, 4, 45], type="l", main="y 5, 45",
         xlab=round(pct.below[45], 3), ylab=round(data$y[4, 45], 2))
    plot(beta.keep[start: i, 1], type="l", main="beta")
    hist.1.range <- range(y.keep[i, , 1], data$y[, 1])
    hist.2.range <- range(y.keep[i, , 7], data$y[, 7])
    hist.3.range <- range(y.keep[i, , 39], data$y[, 39])
    hist.4.range <- range(y.keep[i, , 45], data$y[, 45])
    hist(y.keep[i, , 1], main="imputed day 1", xlim=hist.1.range)
    hist(y.keep[i, , 7], main="imputed day 7", xlim=hist.2.range)
    hist(y.keep[i, , 39], main="imputed day 39", xlim=hist.3.range)
    hist(y.keep[i, , 45], main="imputed day 45", xlim=hist.4.range)
    plot(lambda.2.keep[start:i], type="l", main=lambda.2.t)
    hist(data$y[ , 1], main="actual day 1", xlim=hist.1.range)
    hist(data$y[ , 7], main="actual day 2", xlim=hist.2.range)
    hist(data$y[ , 39], main="actual day 8", xlim=hist.3.range)
    hist(data$y[ , 45], main="actual day 45", xlim=hist.4.range)
    print(paste("iter", i))
  }
}

# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
    vars <- codetools::findGlobals(f)
    found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
    names <- names(found)[found]

    if ((length(names) > 0)) {
      sum.nfncs <- 0
      for (i in 1:length(names)) {
        if(!is.function(eval(parse(text=names[i])))) {sum.nfncs <- sum.nfncs + 1}
      }
      if (sum.nfncs > 0) {
        warning("global variables used: ", paste(names(found)[found], collapse=', '))
        return(invisible(FALSE))
      }
    }

    !any(found)
}

checkStrict(imputeY)
checkStrict(updateBeta)
checkStrict(updateRhoNu)
checkStrict(updateGamma)
checkStrict(updateRhoNuGamma)
checkStrict(updateKnotsTS)
checkStrict(updateLambda1)
checkStrict(updateLambda2)
checkStrict(updatePhiTS)
checkStrict(updateTauGaus)
checkStrict(updateTau)
checkStrict(updateTauTS)
checkStrict(updateZ)
checkStrict(updateZTS)
checkStrict(mcmc)
checkStrict(predictY)
checkStrict(transform$logit)
checkStrict(transform$inv.logit)
checkStrict(transform$probit)
checkStrict(transform$inv.probit)
checkStrict(transform$copula)
checkStrict(transform$inv.copula)
checkStrict(dhn)
checkStrict(phn)
checkStrict(qhn)
checkStrict(mhupdate)
checkStrict(rss)
checkStrict(rTNorm)
checkStrict(CorFx)
checkStrict(eig.inv)
checkStrict(chol.inv)
checkStrict(mem)
checkStrict(makeKnotsTS)
checkStrict(makeTauTS)
checkStrict(makeZTS)
checkStrict(rpotspatTS)
checkStrict(QuantScore)
checkStrict(BrierScore)