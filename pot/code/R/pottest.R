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
nt <- 15
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

beta.t <- c(10, 0, 0)
sigma.t <- 1 / rgamma(nt, 1, 1)
delta.t <- 0.99
rho.t <- 0.1
nu.t <- 0.5
alpha.t <- 0.9

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                   rho=rho.t, nu=nu.t, alpha=alpha.t)

y <- data$y

hist(y)

thresh <- 0.90