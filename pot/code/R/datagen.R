rm(list=ls())

library(fields)
library(geoR)

source("auxfunctions_new.R")
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
tau.alpha.t <- 100
tau.beta.t  <- 100

alpha.t <- 5
rho.t <- 0.1
nu.t  <- 0.5
s2n.t <- 0.8
nknots <- 3

data <- RPotSpatial(nt, s, x, beta=beta.t, tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                    alpha=alpha.t, rho=rho.t, nu=nu.t, s2n=s2n.t, nknots=nknots)
hist(data$y)
data$tau.knots
