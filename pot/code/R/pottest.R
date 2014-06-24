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
nu.t <- 0.75
delta.t <- 0.5
sigma.t <- rep(1, nt)
alpha.t <- 0.9

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=15000, burn=10000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=0.5,
            nu.init=0.5, alpha.init=0.5, delta.init=0,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, fixsigma=T, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=F)
# seems to work alright to get in the neighborhood of rho, nu, alpha, and delta


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

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0
sigma.t <- rep(10, nt)
alpha.t <- 0

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=15000, burn=10000, update=1000, iterplot=T,
            beta.init=c(10, 5, 7), sigma.init=rep(1, nt), rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=F, fixsigma=T, 
            fixrho=T, fixnu=T, fixalpha=T, fixdelta=T)
#estimates beta pretty well

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
nt <- 30
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.9
sigma.t <- 1 / rgamma(nt, 1, 1)
# sigma.t <- rep(1, nt)
alpha.t <- 0.8

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)
# sigma1 = 4.807, sigma 13 = 8.891
# z11 = 0.348, z13 = 2.177

fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=20000, burn=15000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=F, fixbeta=F, fixsigma=F, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=F)

rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

set.seed(1)

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

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.1
sigma.t <- 1 / rgamma(nt, 1, 1)
# sigma.t <- rep(1, nt)
alpha.t <- 0.8

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)
# sigma1 = 4.807, sigma 13 = 8.891
# z11 = 0.348, z13 = 2.177
# z11 = 0.237, z1,16=2.440, z1,21=2.012, z1,30=8.754

fit <- mcmc(y=y, s=s, x=x, thresh=0.9, nknots=1,
            iters=10000, burn=5000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=0.5,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=F, fixbeta=F, fixsigma=F, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=F)
# works reasonably well provided that delta isn't too high.
# if delta is close to 1 and thresh is also close to 1, then there are issues 
# estimating rho.
# for the time being, going to do an unthresholded analysis to estimate rho first.
# then fixing rho at the posterior median.


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
nt <- 30
nsets <- 1
nknots <- 1

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.9
sigma.t <- 1 / rgamma(nt, 1, 1)
# sigma.t <- rep(1, nt)
alpha.t <- 0.8

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=1)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)
# sigma1 = 4.807, sigma 13 = 8.891
# z11 = 0.348, z13 = 2.177

start.time.1 <- proc.time()
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=1,
            iters=20000, burn=15000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=F, fixbeta=F, fixsigma=F, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=F)
end.time.1 <- proc.time()

end.time.1 - start.time.1

#    user  system elapsed 
# 431.873   4.738 434.849 

#### Multiple knots per day

rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

set.seed(1)

# iid n(0, 1)
# data settings
s <- cbind(runif(50), runif(50))
ns <- nrow(s)
nt <- 30
nsets <- 1
nknots <- 3

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.5
sigma.t <- 1 / rgamma(nt, 1, 1)
# sigma.t <- rep(1, nt)
alpha.t <- 0.9

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, sigma=sigma.t, delta=delta.t,
                    rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=nknots)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
# hist(y, breaks=30)
# sigma1 = 1.218, sigma 13 = 10.444
# z1,3 = 6.156, z3,6=0.582, z3,16=15.000, z1,30=4.700

start.time.3 <- proc.time()
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=nknots,
            iters=20000, burn=15000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, knots.init=knots.t, z.init=z.knots.t,
            fixknots=F, fixz=F, fixbeta=F, fixsigma=F, 
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=F)
end.time.3 <- proc.time()

#     user   system  elapsed 
# 1132.833    5.509 1133.624 


#### estimating hyperparams for sigma (one knot per day)

rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

set.seed(1)

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

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.5

sigma.alpha.t <- 5
sigma.beta.t <- 2

alpha.t <- 0.9

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
# sigma5 = 0.382, sigma 30 = 1.218
# z1,3 = 6.156, z3,6=0.582, z3,16=15.000, z1,30=4.700

start.time <- proc.time()
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=nknots,
            iters=10000, burn=5000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.knots.t, 
            sigma.alpha.init=sigma.alpha.t, 
            sigma.beta.init=sigma.beta.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, 
            knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, 
            fixsigma=F, fixsigma.alpha=F, fixsigma.beta=F,
            fixrho=T, fixnu=T, fixalpha=T, fixdelta=T,
            sigma.by.knots=T)
end.time <- proc.time()

end.time - start.time 

# Using MH sampling
#   user  system elapsed 
# 45.747   0.706  50.456 

# Using Gibbs sampling
#   user  system elapsed 
# 31.108   0.479  31.844 

rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

source("auxfunctions.R")
source("mcmc.R")

# set.seed(1)

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

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.5

sigma.alpha.t <- 5
sigma.beta.t <- 2

alpha.t <- 0.9

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
# sigma5 = 0.382, sigma 30 = 1.218
# z1,3 = 6.156, z3,6=0.582, z3,16=15.000, z1,30=4.700

start.time <- proc.time()
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=nknots,
            iters=10000, burn=5000, update=1000, iterplot=T,
            beta.init=beta.t, sigma.init=sigma.knots.t, 
            sigma.alpha.init=sigma.alpha.t, 
            sigma.beta.init=sigma.beta.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            debug=F, 
            knots.init=knots.t, z.init=z.knots.t,
            fixknots=F, fixz=F, fixbeta=F, 
            fixsigma=F, fixsigma.alpha=F, fixsigma.beta=F,
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=F,
            sigma.by.knots=T)
end.time <- proc.time()

end.time - start.time 
# Timing: Nothing fixed, 1 knot, 30 days, 50 sites, 0 threshold, 10000 iters, 5000 burn.
#    user  system elapsed 
# 180.361   1.681 181.302 

rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)

source("auxfunctions.R")
source("mcmc.R")

# set.seed(1000)

# iid n(0, 1)
# data settings
s <- cbind(runif(50), runif(50))
ns <- nrow(s)
nt <- 30
nsets <- 1
nknots <- 3

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0

tau.alpha.t <- 3
tau.beta.t <- 1

alpha.t <- 0.5

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, 
                    tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                    delta=delta.t, rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=nknots)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
tau.knots.t <- data$tau.knots
tau.sites.t <- data$tau.sites
sigma.knots.t <- 1 / tau.knots.t
sigma.sites.t <- 1 / tau.sites.t

# hist(y, breaks=30)

# tau1,1 = 2.923, tau3,20 = 0.510
# z1,8 = 1.208, z3,21 = 0.008

source("auxfunctions.R")
source("mcmc.R")

start.time <- proc.time()
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=nknots,
            iters=15000, burn=10000, update=1000, iterplot=T,
            beta.init=beta.t, tau.init=tau.knots.t, 
            tau.alpha.init=tau.alpha.t, 
            tau.beta.init=tau.beta.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            tau.beta.a=0.01, tau.beta.b=0.01,
            debug=F, 
            knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, 
            fixtau=F, fixtau.alpha=F, fixtau.beta=F,
            fixrho=T, fixnu=T, fixalpha=T, fixdelta=T,
            tau.by.knots=T)
end.time <- proc.time()

# works when alpha = 0

rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)

source("auxfunctions.R")
source("mcmc.R")

# set.seed(1000)

# iid n(0, 1)
# data settings
s <- cbind(runif(50), runif(50))
ns <- nrow(s)
nt <- 90
nsets <- 1
nknots <- 2

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  x[, t, 2] <- s[, 1]
  x[, t, 3] <- s[, 2]
}

beta.t <- c(10, -2, 3)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0

tau.alpha.t <- 3
tau.beta.t <- 1

alpha.t <- 0.8

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, 
                    tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                    delta=delta.t, rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=nknots)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
tau.knots.t <- data$tau.knots
tau.sites.t <- data$tau.sites
sigma.knots.t <- 1 / tau.knots.t
sigma.sites.t <- 1 / tau.sites.t

source("auxfunctions.R")
source("mcmc.R")

start.time <- proc.time()
fit <- mcmc(y=y, s=s, x=x, thresh=0, nknots=nknots,
            iters=10000, burn=5000, update=100, iterplot=T,
            beta.init=beta.t, tau.init=tau.knots.t, 
            tau.alpha.init=tau.alpha.t, 
            tau.beta.init=tau.beta.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            tau.beta.a=0.01, tau.beta.b=0.01,
            debug=F, 
            knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, 
            fixtau=F, fixtau.alpha=F, fixtau.beta=F,
            fixrho=F, fixnu=F, fixalpha=F, fixdelta=T,
            tau.by.knots=T)
end.time <- proc.time()

# hist(y, breaks=30)

# seed 1000 - 30 days
# tau1,1 = 2.923, tau3,20 = 0.510

# seed 1000 - 50 days
# tau1,1 = 4.163, tau3,20 = 2.311

# seed 1001 - 30 days
# tau1,1 = 3.826, tau3,20 = 0.592

# seed 1001 - 50 days
# tau1,1 = 1.148, tau3,20 = 5.609

# seed 1002 - 30 days
# tau1,1 = 4.500, tau3,20 = 1.715

# seed 1002 - 50 days
# tau1,1 = 3.697, tau3,20 = 5.346

# for debugging mcmc
y=y; s=s; x=x; thresh=0; nknots=nknots
thresh.quant=T; scale=T
s.pred=NULL; x.pred=NULL
iters=15000; burn=10000; update=1000; iterplot=F
beta.init=beta.t; tau.init=tau.knots.t
tau.alpha.init=tau.alpha.t
tau.beta.init=tau.beta.t; rho.init=rho.t
nu.init=nu.t; alpha.init=alpha.t; delta.init=delta.t
tau.beta.a=0.01; tau.beta.b=0.01
logrho.m=-2; logrho.s=1
lognu.m=-1.2; lognu.s=1
alpha.m=0; alpha.s=1
debug=F 
knots.init=knots.t; z.init=z.knots.t
fixknots=T; fixz=T; fixbeta=T
fixtau=F; fixtau.alpha=F; fixtau.beta=F
fixrho=F; fixnu=F; fixalpha=F; fixdelta=T
tau.by.knots=T

end.time - start.time 

####  Checking to make sure prediction functions work
rm(list=ls())

library(fields)
library(geoR)
library(mvtnorm)

source("auxfunctions.R")
source("mcmc.R")

set.seed(1000)

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

beta.t <- c(10, 0, 0)
rho.t <- 0.1
nu.t <- 0.5
delta.t <- 0.5

tau.alpha.t <- 3
tau.beta.t <- 1

alpha.t <- 0.8

data <- rpotspatial(nt=nt, s=s, x=x, beta=beta.t, 
                    tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                    delta=delta.t, rho=rho.t, nu=nu.t, alpha=alpha.t, nknots=nknots)

y <- data$y
z.knots.t <- data$z.knots
z.sites.t <- data$z.sites
knots.t <- data$knots
tau.knots.t <- data$tau.knots
tau.sites.t <- data$tau.sites
sigma.knots.t <- 1 / tau.knots.t
sigma.sites.t <- 1 / tau.sites.t

source("auxfunctions.R")
source("mcmc.R")

np <- 10
s.pred <- matrix(runif(np * 2), np, 2)
x.pred <- array(1, c(np, nt, 3))
for (t in 1:nt) {
  x.pred[, t, 2] <- s.pred[, 1]
  x.pred[, t, 3] <- s.pred[, 2]
}

# Rprof(line.profiling=T)
start.time <- proc.time()
fit <- mcmc(y=y, s=s, x=x, # s.pred=s.pred, x.pred=x.pred,
            thresh=0, nknots=nknots,
            iters=5000, burn=1000, update=100, iterplot=T,
            beta.init=beta.t, tau.init=tau.knots.t, 
            tau.alpha.init=tau.alpha.t, 
            tau.beta.init=tau.beta.t, rho.init=rho.t,
            nu.init=nu.t, alpha.init=alpha.t, delta.init=delta.t,
            tau.beta.a=0.01, tau.beta.b=0.01,
            debug=F, 
            knots.init=knots.t, z.init=z.knots.t,
            fixknots=T, fixz=T, fixbeta=T, 
            fixtau=F, fixtau.alpha=F, fixtau.beta=F,
            fixrho=T, fixnu=F, fixalpha=F, fixdelta=T,
            tau.by.knots=T)
end.time <- proc.time()
# Rprof(NULL)
summaryRprof(lines="show")
# for debugging mcmc
y=y; s=s; x=x; s.pred=s.pred; x.pred=x.pred; 
thresh=0; nknots=nknots
thresh.quant=T; scale=T
iters=15000; burn=10000; update=1000; iterplot=F
beta.init=beta.t; tau.init=tau.knots.t
tau.alpha.init=tau.alpha.t
tau.beta.init=tau.beta.t; rho.init=rho.t
nu.init=nu.t; alpha.init=alpha.t; delta.init=delta.t
tau.beta.a=0.01; tau.beta.b=0.01
logrho.m=-2; logrho.s=1
lognu.m=-1.2; lognu.s=1
alpha.m=0; alpha.s=1
debug=F 
knots.init=knots.t; z.init=z.knots.t
fixknots=T; fixz=T; fixbeta=T
fixtau=F; fixtau.alpha=F; fixtau.beta=F
fixrho=F; fixnu=F; fixalpha=F; fixdelta=T
tau.by.knots=T