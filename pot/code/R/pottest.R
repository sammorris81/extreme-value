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

beta.t <- c(0, 0, 0)
sigma.t <- vector(mode="list", length=4)
for (i in 1:4) {
  sigma.t[[i]] <- 1 / rgamma(nt, 4, 4)
  # sigma.t[[i]] <- rep(1, nt)
}

rho.t <- 0.1
nu.t <- 0.5
alpha.t <- 1

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
             fixknots=F, fixz=F, fixbeta=T, fixsigma=F, 
             fixrho=T, fixnu=T, fixalpha=T, fixdelta=T)

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
             fixrho=T, fixnu=T, fixalpha=T, fixdelta=F)
 
 
# testing this out from simple to complicated                   
rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

set.seed(2087)

# iid n(0, 1)
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

beta.t <- c(0, 0, 0)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0
sigma.t <- rep(1, nt)
alpha.t <- 0

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=10000, burn=5000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, fixsigma=F, 
            fixrho=T, fixnu=T, fixalpha=T, fixdelta=T)
# works pretty well


rm(list=ls())
source("auxfunctions.R")
source("mcmc.R")
# iid n(0, sigma2)
# data settings
set.seed(2087)
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

beta.t <- c(0, 0, 0)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0
sigma.t <- 1 / rgamma(nt, 1, 1)
alpha.t <- 0

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots

# sigma1 = 0.8379, sigma3 = 23.7572
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=10000, burn=5000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, fixsigma=F, 
            fixrho=T, fixnu=T, fixalpha=T, fixdelta=T)
# seems to work well now


rm(list=ls())
library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")
# data settings
set.seed(2087)
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

beta.t <- c(0, 0, 0)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0
# sigma.t <- 1 / rgamma(nt, 10, 10)
sigma.t <- rep(5, nt)
alpha.t <- 0.9

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots

# sigma1 = 0.8317, sigma3 = 1.5020
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=20000, burn=15000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, fixsigma=F, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=T)
# works alright

               
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
nt <- 100
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(0, 0, 0)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.7
sigma.t <- rep(10, nt)
alpha.t <- 0.8

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=7000, burn=4000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=0,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, fixsigma=T, 
            fixrho=T, fixnu=T, fixalpha=T, fixdelta=F)
            
# Works for delta alone


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
nt <- 100
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(0, 0, 0)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.7
sigma.t <- rep(10, nt)
alpha.t <- 0.8

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=7000, burn=4000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, fixsigma=T, 
            fixrho=T, fixnu=T, fixalpha=F, fixdelta=F)
# seems to work


rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

# set.seed(1234)

# iid n(0, 1)
# data settings
s <- cbind(runif(50), runif(50))
ns <- nrow(s)
nt <- 100
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(0, 0, 0)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0
sigma.t <- rep(1, nt)
alpha.t <- 0.6

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=15000, burn=10000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=0.5, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, fixsigma=T, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=T)