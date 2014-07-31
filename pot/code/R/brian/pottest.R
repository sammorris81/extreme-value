rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

set.seed(1234)

# iid n(0, 1)
# data settings
s <- cbind(runif(50), runif(50))
ns <- nrow(s)
# nt <- 100  # Working
nt <- 30
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

# beta.t <- c(0, 0, 0)  # Working
beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.9
alpha.t <- 0.8

#sigma.t <- rep(1, nt)
# sigma.t <- rep(10,nt)  # Working
sigma.t <- 1 / rgamma(nt, 1, 1)
data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)
# hist(data$y,breaks=25)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=20000, burn=15000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=F, fixz=F, fixbeta=F, fixsigma=F, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=F)
