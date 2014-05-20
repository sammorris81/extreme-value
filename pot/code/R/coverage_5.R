rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)

source("auxfunctions.R")
source("mcmc.R")

set.seed(1234)

delta.t <- 0.5
image.name <- "coverage_del_50.RData"

# iid n(0, 1)
# data settings
s <- cbind(runif(50), runif(50))
ns <- nrow(s)
nt <- 30
nsets <- 40
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
# sigma.t <- rep(1, nt)
alpha.t <- 0.8

fit.store <- vector(mode="list", length=nsets)
data.store <- vector(mode="list", length=nsets)
sigma.store <- matrix(NA, nrow=nt, ncol=nsets)

coverage.z <- matrix(0, nrow=nknots, ncol=nt)
coverage.sigma <- rep(0, nt)
coverage.beta <- rep(0, 3)
coverage.spatcov <- rep(0, 4) # delta, rho, nu, alpha

iters <- 15000
burn <- 10000
q.start <- burn + 1

for (set in 1:nsets) {
  set.seed(set)
  sigma.t <- 1 / rgamma(nt, 1, 1)
  data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                      rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=nknots)

  y <- data$y
  z.knots.t <- data$z.knots
  z.sites.t <- data$z.sites
  knots.t <- data$knots
# hist(y, breaks=30)
# sigma1 = 4.807, sigma 13 = 8.891
# z11 = 0.348, z13 = 2.177

  fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
              iters=iters, burn=burn, update=1000, iterplot=F, 
              beta.init=c(0, 0, 0), sigma.init=rep(1, nt), rho.init=0.5,
              nu.init=0.5, alpha.init=0.5, delta.init=0,
              debug=F, knots.init=knots.t, z.init=z.knots.t,
              fixknots=T, fixz=F, fixbeta=F, fixsigma=F, 
              fixrho=F, fixnu=F, fixalpha=F, fixdelta=F)
              
  for (k in 1:nknots) { for (t in 1:nt) {
    quants <- quantile(fit$z[q.start:iters, k, t], probs=c(0.025, 0.975))
    coverage.temp <- z.knots.t[k, t] > quants[1] & z.knots.t[k, t] < quants[2]
    coverage.z[k, t] <- coverage.z[k, t] + coverage.temp
  }}
  
  for (t in 1:nt) {
    quants <- quantile(fit$sigma[q.start:iters, t], probs=c(0.025, 0.975))
    coverage.temp <- sigma.t[t] > quants[1] & sigma.t[t] < quants[2]
    coverage.sigma[t] <- coverage.sigma[t] + coverage.temp
  }
  
  for (p in 1:3) {
    quants <- quantile(fit$beta[q.start:iters, p], probs=c(0.025, 0.975))
    coverage.temp <- beta.t[p] > quants[1] & beta.t[p] < quants[2]
    coverage.beta[p] <- coverage.beta[p] + coverage.temp
  }
  
  quants <- quantile(fit$delta[q.start:iters], probs=c(0.025, 0.975))
  coverage.temp <- delta.t > quants[1] & delta.t < quants[2]
  coverage.spatcov[1] <- coverage.spatcov[1] + coverage.temp
  
  quants <- quantile(fit$rho[q.start:iters], probs=c(0.025, 0.975))
  coverage.temp <- rho.t > quants[1] & rho.t < quants[2]
  coverage.spatcov[2] <- coverage.spatcov[2] + coverage.temp
  
  quants <- quantile(fit$nu[q.start:iters], probs=c(0.025, 0.975))
  coverage.temp <- nu.t > quants[1] & nu.t < quants[2]
  coverage.spatcov[3] <- coverage.spatcov[3] + coverage.temp
  
  quants <- quantile(fit$alpha[q.start:iters], probs=c(0.025, 0.975))
  coverage.temp <- alpha.t > quants[1] & alpha.t < quants[2]
  coverage.spatcov[4] <- coverage.spatcov[4] + coverage.temp
  
  data.store[[set]] <- data
  fit.store[[set]]  <- fit
  sigma.store[, set] <- sigma.t
  
  save.image(image.name)
  cat("Dataset", set, "finished \n")
}

coverage.z <- coverage.z / nsets
coverage.sigma <- coverage.sigma / nsets
coverage.beta <- coverage.beta / nsets
coverage.spatcov <- coverage.spatcov / nsets

save.image(image.name)
