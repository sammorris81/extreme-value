rm(list=ls())

options(warn=2)
library(fields)
library(geoR)
library(mvtnorm)

source("auxfunctions.R")
source("mcmc.R")

set.seed(1234)

delta.t <- 0.5
image.name <- "coverage_del_50.RData"

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
alpha.t <- 0.8

sigma.alpha.t <- 5
sigma.beta.t <- 2

fit.store <- vector(mode="list", length=nsets)
data.store <- vector(mode="list", length=nsets)

coverage.z <- matrix(0, nrow=nknots, ncol=nt)
coverage.sigma <- matrix(0, nrow=nknots, ncol=nt)
coverage.sigma.hyper <- rep(0, 2)  # sigma.alpha, sigma.beta
coverage.beta <- rep(0, 3)
coverage.spatcov <- rep(0, 4) # delta, rho, nu, alpha

iters <- 15000
burn <- 10000
q.start <- burn + 1

tic <- proc.time()
for (set in 1:nsets) {
  set.seed(set)
  data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, 
                      sigma.alpha=sigma.alpha.t, sigma.beta=sigma.beta.t, 
                      delta=delta.t, rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=nknots)

  y <- data$y
  z.knots.t <- data$z.knots
  z.sites.t <- data$z.sites
  knots.t <- data$knots
  sigma.knots.t <- data$sigma.knots
  sigma.sites.t <- data$sigma.sites
  
# hist(y, breaks=30)
# sigma1 = 4.807, sigma 13 = 8.891
# z11 = 0.348, z13 = 2.177
  cat("Dataset", set, "started \n")

  fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
              iters=iters, burn=burn, update=1000, iterplot=F, 
              beta.init=c(0, 0, 0), sigma.init=1, 
              sigma.alpha.init=1, sigma.beta.init=1,
              rho.init=0.5, nu.init=0.5, alpha.init=0.5, delta.init=0,
              debug=F, knots.init=knots.t, z.init=z.knots.t,
              fixknots=F, fixz=F, fixbeta=F, 
              fixsigma=F, fixsigma.alpha=F, fixsigma.beta=F,
              fixrho=F, fixnu=F, fixalpha=F, fixdelta=F,
              sigma.by.knots=T)
              
  for (k in 1:nknots) { for (t in 1:nt) {
    quants <- quantile(fit$z[q.start:iters, k, t], probs=c(0.025, 0.975))
    coverage.temp <- z.knots.t[k, t] > quants[1] & z.knots.t[k, t] < quants[2]
    coverage.z[k, t] <- coverage.z[k, t] + coverage.temp
  }}
  
  for (k in 1:nknots) { for (t in 1:nt) {
    quants <- quantile(fit$sigma[q.start:iters, k, t], probs=c(0.025, 0.975))
    coverage.temp <- sigma.knots.t[k, t] > quants[1] & sigma.knots.t[k, t] < quants[2]
    coverage.sigma[t] <- coverage.sigma[k, t] + coverage.temp
  }}
  
  quants <- quantile(fit$sigma.alpha[q.start:iters], probs=c(0.025, 0.975))
  coverage.temp <- sigma.alpha.t > quants[1] & sigma.alpha.t < quants[2]
  coverage.sigma.hyper[1] <- coverage.sigma.hyper[1] + coverage.temp
  
  quants <- quantile(fit$sigma.beta[q.start:iters], probs=c(0.025, 0.975))
  coverage.temp <- sigma.beta.t > quants[1] & sigma.beta.t < quants[2]
  coverage.sigma.hyper[2] <- coverage.sigma.hyper[2] + coverage.temp
  
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
  
  toc <- proc.time()
  
  save.image(image.name)
  cat("Dataset", set, "finished \n")
}

coverage.z <- coverage.z / nsets
coverage.sigma <- coverage.sigma / nsets
coverage.sigma.hyper <- coverage.sigma.hyper / nsets
coverage.beta <- coverage.beta / nsets
coverage.spatcov <- coverage.spatcov / nsets

save.image(image.name)
