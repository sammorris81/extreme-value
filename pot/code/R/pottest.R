options(warn=2)
#set.seed(0820)
source("mcmc.R")
source("auxfunctions.R")

ns <- 20
nt <- 50

s <- cbind(seq(0, 1, length=ns), 0.5)

a <- 2
b <- 8

rho    <- .1
nu     <- 0.5
alpha  <- 0.95
mu     <- 15    
d <- as.matrix(dist(s))
C <- CorFx(d=d, alpha=alpha, rho=rho, nu=nu)

nk <- 5
g <- rep(1:4, each=5)
tau <- matrix(rgamma(nk * nt, a, b), nk, nt)
taug <- tau[g, ]

nk <- 1
g <- rep(1, 20)
tau <- matrix(rgamma(nt, a, b), nk, nt)
taug <- tau[g, ]

nk <- 1
g <- rep(1, 20)
tau <- matrix(rgamma(1, a, b), nk, nt)
taug <- tau[g, ]


sd <- 1 / sqrt(taug)
Y <- t(chol(C)) %*% matrix(rnorm(ns * nt), ns, nt)
Y <- mu + sd * Y


x <- array(1,c(ns,nt,1))

sp <- s
sp[, 2] <- 0.6
xp <- x
dp <- as.matrix(dist(sp))
Cp <- CorFx(d=d, alpha=alpha, rho=rho, nu=nu) 
yp <- t(chol(Cp)) %*% matrix(rnorm(ns * nt), ns, nt)
yp <- mu + sd * yp

thresh <- 0

fit<-mcmc(Y, s, x, s.pred=sp, x.pred=xp, method="t",
          thresh=thresh, thresh.quant=T, nknots=nk, iterplot=T,
          iters=5000, burn=1000, update=100, thin=1,
          rho.init=rho, nu.init=nu, alpha.init=alpha)

# plot posterior predictions
lower <- apply(fit$yp, c(2, 3), quantile, probs=0.025)
upper <- apply(fit$yp, c(2, 3), quantile, probs=0.975)

days <- 16
par(mfrow=c(4, 4))
for (day in 1:days) {
  plot(yp[, day], pch=19, ylim=c(min(lower[, day]), max(upper[, day])), main=paste("day ", day, sep=""))
  lines(lower[, day], lty=2)
  lines(upper[, day], lty=2)
}
for (day in 17:(2*days)) {
  plot(yp[, day], pch=19, ylim=c(min(lower[, day]), max(upper[, day])), main=paste("day ", day, sep=""))
  lines(lower[, day], lty=2)
  lines(upper[, day], lty=2)
}
for (day in 33:(3*days)) {
  plot(yp[, day], pch=19, ylim=c(min(lower[, day]), max(upper[, day])), main=paste("day ", day, sep=""))
  lines(lower[, day], lty=2)
  lines(upper[, day], lty=2)
}

# Add in skew
options(warn=2)
rm(list=ls())
source("mcmc.R")
source("auxfunctions.R")

ns <- 20
nt <- 50

s <- cbind(seq(0, 1, length=ns), 0.5)

a <- 2
b <- 8

rho    <- 0.1
nu     <- 0.5
alpha  <- 0.95
x.beta <- 15    
d <- as.matrix(dist(s))
C <- CorFx(d=d, alpha=alpha, rho=rho, nu=nu)

nk <- 4
g <- rep(1:4, each=5)
tau <- matrix(rgamma(nk * nt, a, b), nk, nt)
taug <- tau[g, ]
z.alpha <- 0

nk <- 1
g <- rep(1, ns)
tau <- matrix(rgamma(nt, a, b), nk, nt)
taug <- tau[g, ]
z.alpha <- 0

nk <- 1
g <- rep(1, ns)
tau <- matrix(rgamma(1, a, b), nk, nt)
taug <- tau[g, ]
z.alpha <- 5


sdg <- 1 / sqrt(taug)
sd <- 1 / sqrt(tau)
Y <- t(chol(C)) %*% matrix(rnorm(ns * nt), ns, nt)
z <- matrix(abs(rnorm(nk * nt, 0, 1)), nk, nt) * sd 
zg <- z[g, ]
mu <- x.beta + z.alpha * zg
Y <- mu + sdg * Y
hist(Y)

x <- array(1,c(ns,nt,1))

sp <- s
sp[, 2] <- 0.6
xp <- x
dp <- as.matrix(dist(sp))
Cp <- CorFx(d=d, alpha=alpha, rho=rho, nu=nu) 
yp <- t(chol(Cp)) %*% matrix(rnorm(ns * nt), ns, nt)
yp <- mu + sdg * yp

thresh <- 0.90

fit<-mcmc(Y, s, x, s.pred=sp, x.pred=xp, method="t", skew=F,
          thresh=0.9, thresh.quant=T, nknots=nk, iterplot=T,
          iters=10000, burn=5000, update=100, thin=1,
          rho.init=rho, nu.init=nu, alpha.init=alpha, z.init=z,
          z.alpha.m=1, z.alpha.s=1, z.alpha.init=0)
          
test <- rnorm(10000, 0, 1)
norm.test <- pnorm(test)
hist(norm.test, breaks=80)


Y <- rpotspat(ns, nt, X=x, s, beta=15, alpha=0.95, gau.rho=0.1, t.rho=0.1, nu=0.5, mixprob=1, z.alpha=3, tau.alpha=2, tau.beta=8, nknots=5)
hist(Y$y)