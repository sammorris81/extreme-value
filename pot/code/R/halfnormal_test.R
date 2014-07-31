rm(list=ls())
eta <- 3
xi <- 1

x <- rnorm(1, 0, 1)
u <- xi + sqrt(eta) * abs(x)
theta.true <- 1 / eta

n <- 100
y <- rnorm(n, u, 1)
sigma2 <- 1
tau <- 1 / sigma2

iters <- 20000
u.reps <- rep(NA, iters)

# Gibbs to see if this is properly specified
for(i in 1:iters){
	xistar <- tau * sum(y) / (theta.true + n * tau)
	thetastar <- theta.true + n * tau
	etastar <- 1 / sqrt(thetastar)
	u.reps[i] <- xistar + etastar * rnorm(1, 0, 1)
}

