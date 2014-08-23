#########################################################################
# A simulation study to determine if a thresholded or skew methods improve
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#   1 - Gaussian
#   2 - t-1
#   3 - t-5
#   4 - skew t-1 (alpha = 3)
#   5 - skew t-5 w/partition (alpha = 3)
#   6 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
#
# analysis methods:
#  1 - Gaussian
#  2 - skew t-1
#  3 - skew t-1 (T = 0.90)
#  4 - skew t-5
#  5 - skew t-5 (T = 0.90)
#	
#########################################################################

# clean out previous variables
rm(list=ls())

# libraries
library(fields)
library(geoR)

# necessary functions
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

# data settings
beta.t <- c(10, 0, 0)
nu.t <- 0.5
alpha.t <- 0.9
mixprob.t <- c(0, 1, 1, 1, 1, 0.5)  # 0: Gaussian, 1: t
nknots.t <- c(1, 1, 5, 1, 5, 1)
gau.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10)
t.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.40)
z.alpha.t <- c(0, 0, 0, 3, 3, 0)
tau.alpha.t <- 2
tau.beta.t  <- 8


# covariate data
s         <- cbind(runif(130), runif(130))
ns        <- nrow(s)  
nt        <- 50
nsets     <- 50
nsettings <- length(mixprob.t) 
ntest     <- 30

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

# Storage for datasets
y <- array(NA, dim=c(ns, nt, nsets, nsettings))
tau.t <- z.t <- knots.t <- vector("list", length=nsettings)

for (setting in 1:nsettings) {
  nknots <- nknots.t[setting]
  tau.t.setting <- array(NA, dim=c(nknots, nt, nsets))
  z.t.setting <- array(NA, dim=c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim=c(nknots, 2, nt, nsets))
  for (set in 1:nsets) {
    set.seed(setting*100 + set)
    data <- rpotspat(nt=nt, x=x, s=s, beta=beta.t, alpha=alpha.t, nu=nu.t,
                     gau.rho=gau.rho.t[setting], t.rho=t.rho.t[setting],
                     mixprob=mixprob.t[setting], z.alpha=z.alpha.t[setting],
                     tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                     nknots=nknots.t[setting])
    
    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }
    
  tau.t[[setting]]   <- tau.t.setting
  z.t[[setting]]     <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
}

# plot the data
s1 <- s2 <- seq(0, 1, length=50)
s <- expand.grid(s1, s2)
X <- array(1, dim=c(nrow(s), 2, 3))
X.temp <- cbind(rep(1, nrow(s)), s[, 1], s[, 2])
X[, 1, ] <- X.temp
X[, 2, ] <- X.temp
y.gau <- rpotspat(nt=2, x=X, s=s, beta=c(10, 0, 0), alpha=alpha.t, nu=nu.t,
                  gau.rho=0.1, t.rho=0.1, mixprob=0, z.alpha=0, tau.alpha=tau.alpha.t,
                  tau.beta=tau.beta.t, nknots=1)
par(mfrow=c(1, 2))
quilt.plot(s[, 1], s[, 2], z=y.gau$y[, 1], nx=50, ny=50, main="Gaussian")
hist(y.gau$y[, 1], main="Histogram", xlab="")

y.t1 <- rpotspat(nt=2, x=X, s=s, beta=c(10, 0, 0), alpha=alpha.t, nu=nu.t,
                  gau.rho=0.1, t.rho=0.1, mixprob=1, z.alpha=0, tau.alpha=tau.alpha.t,
                  tau.beta=tau.beta.t, nknots=1)
par(mfrow=c(1, 2))
quilt.plot(s[, 1], s[, 2], z=y.t1$y[, 1], nx=50, ny=50, main="t, K=1")
hist(y.t1$y[, 1], main="Histogram", xlab="")

y.t3 <- rpotspat(nt=2, x=X, s=s, beta=c(10, 0, 0), alpha=alpha.t, nu=nu.t,
                  gau.rho=0.1, t.rho=0.1, mixprob=1, z.alpha=0, tau.alpha=tau.alpha.t,
                  tau.beta=tau.beta.t, nknots=3)
par(mfrow=c(1, 2))
quilt.plot(s[, 1], s[, 2], z=y.t3$y[, 1], nx=50, ny=50, main="t, K=3")
hist(y.t3$y[, 1], main="Histogram", xlab="")

y.st1 <- rpotspat(nt=2, x=X, s=s, beta=c(10, 0, 0), alpha=alpha.t, nu=nu.t,
                  gau.rho=0.1, t.rho=0.1, mixprob=1, z.alpha=3, tau.alpha=tau.alpha.t,
                  tau.beta=tau.beta.t, nknots=1)
par(mfrow=c(1, 2))
quilt.plot(s[, 1], s[, 2], z=y.st1$y[, 1], nx=50, ny=50, main="skew-t, K=1, alpha=3")
hist(y.st1$y[, 1], main="Histogram", xlab="")

y.st3 <- rpotspat(nt=2, x=X, s=s, beta=c(10, 0, 0), alpha=alpha.t, nu=nu.t,
                  gau.rho=0.1, t.rho=0.1, mixprob=1, z.alpha=1, tau.alpha=tau.alpha.t,
                  tau.beta=tau.beta.t, nknots=3)
par(mfrow=c(1, 2))
quilt.plot(s[, 1], s[, 2], z=y.s35$y[, 1], nx=50, ny=50, main="skew-t, K=3, alpha=1")
hist(y.st5$y[, 1], main="Histogram", xlab="")

par(mfrow=c(1, 2))
quilt.plot(s[, 1], s[, 2], z=y.gau$y[, 1], nx=50, ny=50, main="Gaussian")
quilt.plot(s[, 1], s[, 2], z=y.st3$y[, 1], nx=50, ny=50, main="skew-t, K=3, alpha=1")

par(mfrow=c(1, 2))
quilt.plot(s[, 1], s[, 2], z=y.gau$y[, 1], nx=50, ny=50, main="Gaussian")
quilt.plot(s[, 1], s[, 2], z=y.t3$y[, 1], nx=50, ny=50, main="t, K=3")

# remove simstudy "truth" settings to help diagnose errors.
save(y, tau.t, z.t, knots.t, ns, nt, s, nsets, 
     x, ntest,  # covariate data that should be the same for all datasets
     file='simdata.RData')
rm(list=ls())
