rm(list=ls())
source("../R/mcmc.R")
source("../R/auxfunctions.R")

# lowerbound for integrated quantile and brier scores
U <- 0.95
quant.grid <- seq(U, 1, 0.001)
thresh.probs <- seq(0.85, U, 0.01)
thresh.true <- 0.90

# settings for GPD and probability of exceedence
set.seed(2087)
n <- 500
nsets <- 15
settings <- 3
X <- matrix(rep(1, n), nrow=n, ncol=1)
# X <- matrix(rnorm(5 * n), n, 5)
# X[, 1] <- 1
beta1 <- c(1, rnorm(4, 0, 1) * rbinom(4, 1, 0.6))
beta2 <- c(1, rnorm(4, 0, 1) * rbinom(4, 1, 0.6))
beta3 <- c(1, rnorm(4, 0, 1) * rbinom(4, 1, 0.6))

# more complicated than we want to see if this is working
# prob.exceed <- expit(X %*% beta1)
# gpd.sig <- exp(as.vector(X %*% beta2))
# gpd.xi <- as.vector(X %*% beta3)

prob.below <- thresh.true
prob.above <- 1 - prob.below
gpd.sig <- 2
gpd.xi <- 0.2

# height for continuity of density
height <- prob.above / gpd.sig

# setting 1: uniform below
lower <- - 1 / height

# setting 2: normal w/sd=1
sig.1 <- 1
mu.1 <- -sqrt(-2 * sig.1^2 * (log(height) + log(sig.1) + 0.5 * log(2 * pi)))
rescale.1 <- pnorm(0, mu.1, sig.1)

# setting 3: normal w/sd=5
sig.2 <- 5
mu.2 <- -sqrt(-2 * sig.2^2 * (log(height) + log(sig.2) + 0.5 * log(2 * pi)))
rescale.2 <- pnorm(0, mu.2, sig.2)

y <- array(NA, dim=c(n, nsets, settings))

# generate data
for (set in 1:nsets) {
  
  for (i in 1:n) {
    u <- runif(1, 0, 1)
    if (u <= prob.below) {
    	y[i, set, 1] <- qunif(u / prob.below, lower, 0)
    } else {
    	# this is different from what happens in POT package
    	adj <- (u - prob.below) / (prob.above)
    	y[i, set, 1] <- qGP(p=adj, scale=gpd.sig, shape=gpd.xi)
    }
  }
  
  for (i in 1:n) {
    u <- runif(1, 0, 1)
    if (u <= prob.below){
    	rescale.u <- u * rescale.1 / prob.below
    	y[i, set, 2] <- qnorm(rescale.u, mu.1, sig.1)
    } else {
    	adj <- (u - prob.below) / (prob.above)
    	y[i, set, 2] <- qGP(p=adj, scale=gpd.sig, shape=gpd.xi)
    }
  }
    
  for (i in 1:n) {
    u <- runif(1, 0, 1)
    if (u <= prob.below){
    	rescale.u <- u * rescale.2 / prob.below
    	y.temp <- qnorm(rescale.u, mu.2, sig.2)
    	if (y.temp > 0) {print("y.temp too big")}
    } else {
    	adj <- (u - prob.below) / (prob.above)
    	y.temp <- qGP(p=adj, scale=gpd.sig, shape=gpd.xi)
    }
    y[i, set, 3] <- y.temp
  }
}


save(U, quant.grid, thresh.probs, thresh.true,
     beta1, beta2, beta3, X,
     gpd.sig, gpd.xi, y,
     file='simdata.RData')
     
rm(list=ls())

# reload the functions for the simulation into the workspace
load(file='simdata.Rdata')
source('../R/mcmc.R')
source('../R/auxfunctions.R')

save.image(file='simdata.RData')


# # ########################################################
# ########    EXAMPLE
# ########################################################
# n <- 500
# X <- matrix(rnorm(5 * n), n, 5)
# X[, 1] <- 1
# prob <- ifelse (X[, 2] > 0, 0.4, 0.1)
# exceed <- exp(rnorm(n, X[, 3], 1))
# y <- ifelse (runif(n) < prob, exceed, 0) + 10
 
fit <- Bayes_GPD(y, X, X, X, thresh=10, iterplot=T, debug=F,
				iters=5000, burn=1000, update=1000)

y.test <- y[, 1, 3]
thresh <- 0
X <- matrix(rep(1, length(y.test)), nrow=length(y.test), ncol=1)
#X[, 2] <- rnorm(length(y), 0, 1)

fit <- Bayes_GPD(y=y.test, X.prob=X, X.sig=X, X.xi=X, 
                 # Xp.prob=X, Xp.sig=X, Xp.xi=X,
                 prob.m=0, prob.s=10, fix.prob=F,
                 sig.m=-2, sig.s=1, fix.sig=F,
                 xi.m=0, xi.s=0.2, fix.xi=F,
                 prob.init=0.5, sig.init=1, xi.init=0.1,
                 thresh=thresh, iters=10000, burn=5000, update=1000,
                 iterplot=T, debug=F)

# # # troubleshooting
# # y <- y - 10
# # y <- y[y > 0]
# # scale <- rep(39, length(y))
# # shape <- rep(0.1, length(y))
# # GP(y, scale, shape)

# # plot densities
# library(POT)
# par(mfrow=c(3, 2))
# xplot <- seq(lower, 10, 0.001)
# yplot <- dunif(xplot, lower, 0) + 
         # dgpd(xplot, 0, gpd.sig, gpd.xi)
# plot(xplot, yplot, type="l", main="uniform below")

# hist(y[, , 1], freq=F, xlim=c(0, max(y[, , 1])))
# xplot <- seq(0, 60, 0.001)
# yplot <- dgpd(xplot, 0, gpd.sig*prob.above, gpd.xi)
# lines(xplot, yplot)

# xplot <- seq(-10, 10, 0.001)
# yplot <- dnorm(xplot, mu.1, sig.1) * (xplot <= 0) + 
         # (prob.above) * dgpd(xplot, 0, gpd.sig, gpd.xi)
# plot(xplot, yplot, type="l", main="normal below with sd=1")

# hist(y[, , 2], freq=F, xlim=c(0, max(y[, , 2])))
# xplot <- seq(0, 10, 0.001)
# yplot <- dgpd(xplot, 0, gpd.sig*prob.above, gpd.xi)
# lines(xplot, yplot)

# xplot <- seq(-10, 10, 0.001)
# yplot <- dnorm(xplot, mu.2, sig.2) * (xplot <= 0) + 
         # (prob.above) * dgpd(xplot, 0, gpd.sig, gpd.xi)
# plot(xplot, yplot, type="l", main="normal below with sd=5")

# hist(y[, , 3], freq=F, xlim=c(0, max(y[, , 3])))
# xplot <- seq(0, 10, 0.001)
# yplot <- dgpd(xplot, 0, gpd.sig*prob.above, gpd.xi)
# lines(xplot, yplot)

# hist(y.1, breaks=100)
# hist(y.2, breaks=100)
# hist(y.3, breaks=100)

# par(mfrow = c(2, 2))
# plot(X[, 2], y)
# boxplot(fit$beta.prob, main="betas in ex prob", outline=F)
# abline(0, 0)
# boxplot(fit$beta.sig, main="betas in GDP scale", outline=F)
# abline(0, 0)
# boxplot(fit$beta.xi, main="betas in GPD shape", outline=F)  
# abline(0, 0)
