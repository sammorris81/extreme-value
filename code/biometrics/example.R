rm(list=ls())

# libraries
library(fields)
library(SpatialTools)

# necessary functions
source('mcmc.R')
source('auxfunctions.R')

set.seed(100)
ns <- 150
nt <- 50
s <- cbind(runif(ns), runif(ns))
X <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
  X[, t, 2] <- s[, 1]
  X[, t, 3] <- s[, 2]
}

# generate from skew-t with 1 knot, lambda = 3, df = 8, no time series
data1 <- rpotspatTS(nt = nt, x = X, s = s, beta = c(10, 0, 0),
                    gamma = 0.9, nu = 0.5, rho = 0.1,
                    phi.z = 0, phi.w = 0, phi.tau = 0,  # time series parameters
                    lambda = 3, tau.alpha = 6, tau.beta = 16,
                    nknots = 1, dist = "t")$y

# fit the model using no time series with 1 knot, skew, and matern covariance.
# Prior distributions:
# beta      ~iid N(0, 20)
# tau.alpha ~ discrete(seq(0.2, 20, by = 0.2))
# tau.beta  ~ gamma(0.1, 0.1)
# rho       ~ uniform(0, 1.2)
# lambda    ~ N(0, 20)
fit <- mcmc(y = data1, s = s, x = X, min.s = c(0, 0), max.s = c(1, 1),
            thresh.all = 0.5, thresh.quant = TRUE, nknots = 1, skew = TRUE,
            method = "t", cov.model = "matern",
            temporalw = FALSE, temporaltau = FALSE, temporalz = FALSE,
            rho.init = 0.5,
            beta.m = 0, beta.s = 20,
            tau.alpha.min = 0.2, tau.alpha.max = 20, tau.alpha.by = 0.2,
            tau.beta.a = 0.1, tau.beta.b = 0.1, rho.upper = 1.2,
            lambda.m = 0, lambda.s = 20,
            iters = 10000, burn = 5000, update = 100, iterplot = TRUE)

# generate from spatiotemporal skew-t model
data2 <- rpotspatTS(nt = nt, x = x, s = s, beta = c(10, 0, 0),
                    gamma = 0.9, nu = 0.5, rho = 0.1,
                    phi.z = 0.5, phi.w = 0.5, phi.tau = 0.5,
                    lambda = 3, tau.alpha = 8, tau.beta = 10,
                    nknots = 1, dist = "t")$y

# fit the model using time series with 1 knot, skew, and matern covariance.
# Prior distributions:
# beta      ~iid N(0, 20)
# tau.alpha ~ discrete(seq(0.2, 20, by = 0.2))
# tau.beta  ~ gamma(0.1, 0.1)
# rho       ~ uniform(0, 1.2)
# lambda    ~ N(0, 20)
fit <- mcmc(y = data1, s = s, x = X, min.s = c(0, 0), max.s = c(1, 1),
            thresh.all = 0.5, thresh.quant = TRUE, nknots = 1, skew = TRUE,
            method = "t", cov.model = "matern",
            temporalw = TRUE, temporaltau = TRUE, temporalz = TRUE,
            rho.init = 0.5,
            beta.m = 0, beta.s = 20,
            tau.alpha.min = 0.2, tau.alpha.max = 20, tau.alpha.by = 0.2,
            tau.beta.a = 0.1, tau.beta.b = 0.1, rho.upper = 1.2,
            lambda.m = 0, lambda.s = 20,
            iters = 10000, burn = 5000, update = 100, iterplot = TRUE)
