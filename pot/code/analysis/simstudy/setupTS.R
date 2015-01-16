# Want to generate a few time series datasets to make sure the results from
# the MCMC are reasonable

rm(list=ls())
library(fields)
library(SpatialTools)
source('../../R/mcmc.R', chdir=T)
source('../../R/auxfunctions.R')

# parameter settings
beta.t      <- c(10, 0, 0)
nu.t        <- 0.5
gamma.t     <- 0.9
nknots.t    <- c(1, 5, 1, 5)
rho.t       <- 1
lambda.t    <- c(0, 0, 3, 3)
tau.alpha.t <- 3
tau.beta.t  <- 8
phi.w.t     <- 0.8
phi.z.t     <- 0.8
phi.tau.t   <- 0.8

# covariates and data setings
ns        <- 144
nt        <- 50
nsets     <- 5
ntest     <- 44
nsettings <- 4

set.seed(20)
s <- cbind(runif(ns, 0, 10), runif(ns, 0, 10))
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
  tau.t.setting   <- array(NA, dim=c(nknots, nt, nsets))
  z.t.setting     <- array(NA, dim=c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim=c(nknots, 2, nt, nsets))
  for (set in 1:nsets) {
    set.seed(setting * 100 + set)
    data <- rpotspatTS(nt=nt, x=x, s=s, beta=beta.t,
                       gamma=gamma.t, nu=nu.t, rho=rho.t,
                       phi.z=phi.z.t, phi.w=phi.w.t, phi.tau=phi.tau.t,
                       lambda=lambda.t[setting], nknots=nknots,
                       tau.alpha=tau.alpha.t, tau.beta=tau.beta.t)
    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }
  tau.t[[setting]]   <- tau.t.setting
  z.t[[setting]]     <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
}

save(y, tau.t, z.t, knots.t, ns, nt, s, nsets, ntest, x, file='simdataTS.RData')