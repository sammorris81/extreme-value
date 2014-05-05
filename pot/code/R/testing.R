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


################################################################
LLike2 <- function(y, x.beta, prec, log.det, z.sites, log=TRUE){
  
  if (missing(y)) {
    stop("y must be defined")
  } 
  ns <- nrow(y)
  nt <- ncol(y)
  if (missing(x.beta)) {
    stop("x.beta must be defined")
  } else if (nrow(x.beta) != ns || ncol(x.beta) != nt) {
    stop("x.beta and y must have conforming size")
  }

  
  if (missing(log.det)) {
    stop("log.det must be defined")
  }
      
  log.like  <- rep(NA, nt)
  
  for (t in 1:nt) {
    y.t <- y[, t]
    mu.t <- x.beta[, t] + z.sites[, t]
    ss.t <- t(y.t - mu.t) %*% prec %*% (y.t - mu.t)
    log.like[t] <- 0.5 * log.det - 0.5 * ss.t
  }
    
  if(!log){
    log.like <- exp(log.like)
  }
    
  return(log.like)
}

y.test <- matrix(y[[4]][,3], 50, 1)
sigma.test <- sigma.t[[4]][1]
delta.test <- delta.t[[4]]
z.sites.test <- matrix(z.sites.t[[4]][, 1], 50, 1)
z.sites.test.sig <- z.sites.test * delta.test
x.beta.test <- matrix(rep(10, 50), 50, 1)
mu.test <- x.beta.test + z.sites.test * delta.test

alpha.t <- .9
nug.sig <- (1 - alpha.t) * (sigma.test * (1 - delta.test^2))
p.sill.t.sig <- alpha.t * (sigma.test * (1 - delta.test^2))

SIG.sig <- varcov.spatial(coords=s, cov.model="matern", nugget=nug.sig, cov.pars=c(p.sill.t.sig, rho.t), kappa=nu.t)$varcov
PREC.sig <- chol2inv(chol(SIG.sig))
log.det.test.sig <- log(det(PREC.sig))


# # dec.sig <- chol(SIG.sig)
# tmp.sig <- forwardsolve(t(dec.sig), y.test - mu.test)
# rss <- colSums(tmp.sig^2)

LLike2(y=y.test, x.beta=x.beta.test, prec=PREC.sig, log.det=log.det.test.sig, z.sites=z.sites.test.sig)



dmvnorm(t(y.test), mu.test, SIG.sig, log=T) + 25 * log(2 * pi)

nug <- (1 - alpha.t)
p.sill.t <- alpha.t
SIG <- varcov.spatial(coords=s, cov.model="matern", nugget=nug, cov.pars=c(p.sill.t, rho.t), kappa=nu.t)$varcov
PREC <- chol2inv(chol(SIG))
log.det.t <- log(det(PREC))

LLike(y=y.test, x.beta=x.beta.test, sigma=sigma.test, delta=delta.test, prec=PREC, log.det=log.det.t, z.sites=z.sites.test, log=T)



