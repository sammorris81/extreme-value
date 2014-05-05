rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

set.seed(2087)

# data settings
s <- cbind(runif(50), runif(50))
ns <- nrow(s)
nt <- 30
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

beta.t <- c(10, 0, 0)
sigma.t <- vector(mode="list", length=4)
for (i in 1:4) {
  sigma.t[[i]] <- 1 / rgamma(nt, 4, 4)
  # sigma.t[[i]] <- rep(1, nt)
}

rho.t <- 0.1
nu.t <- 0.5
alpha.t <- 0

# making sure the data generated looks reasonable
y         <- vector(mode="list", length=4)
z.knots.t <- vector(mode="list", length=4)
z.sites.t <- vector(mode="list", length=4)
knots.t   <- vector(mode="list", length=4)
delta.t   <- c(0, 0.5, 0.9, 0.95)

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t[[1]], delta=delta.t[1],
                   rho=rho.t, nu=nu.t, alpha=alpha.t)          
y[[1]] <- data$y
z.knots.t[[1]] <- data$z.knots
z.sites.t[[1]] <- data$z.sites
knots.t[[1]] <- data$knots
# hist(y[[1]], main="delta = 0")

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t[[2]], delta=delta.t[2],
                   rho=rho.t, nu=nu.t, alpha=alpha.t)          
y[[2]] <- data$y
z.knots.t[[2]] <- data$z.knots
z.sites.t[[2]] <- data$z.sites
knots.t[[2]] <- data$knots
# hist(y[[2]], main="delta = 0.5")

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t[[3]], delta=delta.t[3],
                   rho=rho.t, nu=nu.t, alpha=alpha.t)
y[[3]] <- data$y
z.knots.t[[3]] <- data$z.knots
z.sites.t[[3]] <- data$z.sites
knots.t[[3]] <- data$knots
# hist(y[[3]], main="delta = 0.9")

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t[[4]], delta=delta.t[4],
                   rho=rho.t, nu=nu.t, alpha=alpha.t)
y[[4]] <- data$y
z.knots.t[[4]] <- data$z.knots
z.sites.t[[4]] <- data$z.sites
knots.t[[4]] <- data$knots
# hist(y[[4]], main="delta = 0.95")

source("auxfunctions.R")
source("mcmc.R")

# testing out the MCMC
# sigma1 = 0.7858, sigma3 = 2.1346
# z11 = 0.3614, z13 = 0.6151
fit1 <- mcmc(y=y[[1]], s=s, x=x, thresh=0, nknots=1,
             iters=10000, burn=5000, update=1000, iterplot=T,
             beta.init=beta.t, sigma.init=sigma.t[[1]], rho.init=rho.t,
             nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t[1],
             debug=F, knots.init=knots.t[[1]], z.init=z.knots.t[[1]],
             fixknots=F, fixz=F, fixbeta=F, fixsigma=F, 
             fixrho=F, fixnu=F, fixalpha=F, fixdelta=T)

# sigma1 = 0.9109, sigma3 = 1.2974
# z11 = 0.8964, z13 = 0.8644
fit2 <- mcmc(y=y[[2]], s=s, x=x, thresh=0, nknots=1,
             iters=10000, burn=5000, update=1000, iterplot=T,
             beta.init=beta.t, sigma.init=sigma.t[[2]], rho.init=rho.t,
             nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t[2],
             debug=F, knots.init=knots.t[[2]], z.init=z.knots.t[[2]],
             fixknots=T, fixz=T, fixbeta=T, fixsigma=T, 
             fixrho=T, fixnu=T, fixalpha=T, fixdelta=F)

# sigma1 = 1.4812, sigma3 = 0.7718
# z11 = 1.1067, z13 = 0.3671
fit3 <- mcmc(y=y[[3]], s=s, x=x, thresh=0, nknots=1,
             iters=10000, burn=5000, update=1000, iterplot=T,
             beta.init=beta.t, sigma.init=sigma.t[[3]], rho.init=rho.t,
             nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t[3],
             debug=F, knots.init=knots.t[[3]], z.init=z.knots.t[[3]],
             fixknots=T, fixz=T, fixbeta=T, fixsigma=T, 
             fixrho=T, fixnu=T, fixalpha=T, fixdelta=F)

# sigma1 = 2.2898, sigma3 = 0.5745
# z11 = 1.4383, z13 = 0.5165
fit4 <- mcmc(y=y[[4]], s=s, x=x, thresh=0, nknots=1,
             iters=10000, burn=5000, update=1000, iterplot=T,
             beta.init=beta.t, sigma.init=sigma.t[[4]], rho.init=rho.t,
             nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t[4],
             debug=F, knots.init=knots.t[[4]], z.init=z.knots.t[[4]],
             fixknots=T, fixz=T, fixbeta=T, fixsigma=T, 
             fixrho=T, fixnu=T, fixalpha=T, fixdelta=T)
 