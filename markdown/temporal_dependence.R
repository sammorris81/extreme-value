rm(list = ls())

library(fields)
library(SpatialTools)
library(evd)

# library(compiler)
# enableJIT(3)

#### Load simdata
source('../code/R/mcmc_cont_lambda.R', chdir=T)
source('../code/R/auxfunctions.R')

options(warn=2)


set.seed(100)
ns <- 5
nt <- 100000
X <- array(1, dim = c(ns, nt, 1))
S <- matrix(runif(ns * 2), ns, 2)
data <- rpotspatTS_cont_lambda(nt = nt, x = X, s = S, beta = 0, gamma = 0.90, 
                               nu = 0.5, rho = 0.2, phi.z = 0.9, phi.w = 0.9, 
                               phi.tau = 0.9, lambda = 3, tau.alpha = 1, 
                               tau.beta = 1, nknots = 1, dist = "t")

orig.s1 <- data$y[1, 1:(nt - 1)]
lag.s1  <- data$y[1, 2:nt]
data.s1 <- cbind(orig.s1, lag.s1)

orig.s2 <- data$y[2, 1:(nt - 1)]
lag.s2  <- data$y[2, 2:nt]
data.s2 <- cbind(orig.s2, lag.s2)

orig.s3 <- data$y[3, 1:(nt - 1)]
lag.s3  <- data$y[3, 2:nt]
data.s3 <- cbind(orig.s3, lag.s3)

orig.s4 <- data$y[4, 1:(nt - 1)]
lag.s4  <- data$y[4, 2:nt]
data.s4 <- cbind(orig.s4, lag.s4)

orig.s5 <- data$y[5, 1:(nt - 1)]
lag.s5  <- data$y[5, 2:nt]
data.s5 <- cbind(orig.s5, lag.s5)

chiplot(data = data.s1, which = 1, main1 = "site 1")
chiplot(data = data.s2, which = 1, main1 = "site 2")
chiplot(data = data.s3, which = 1, main1 = "site 3")
chiplot(data = data.s4, which = 1, main1 = "site 4")
chiplot(data = data.s5, which = 1, main1 = "site 5")

chiplot(data = data.s1, which = 2, main2 = "site 1")
chiplot(data = data.s2, which = 2, main2 = "site 2")
chiplot(data = data.s3, which = 2, main2 = "site 3")
chiplot(data = data.s4, which = 2, main2 = "site 4")
chiplot(data = data.s5, which = 2, main2 = "site 5")

set.seed(100)
ns <- 5
nt <- 100000
X <- array(1, dim = c(ns, nt, 1))
S <- matrix(runif(ns * 2), ns, 2)
data <- rpotspatTS_cont_lambda(nt = nt, x = X, s = S, beta = 0, gamma = 0.90, 
                               nu = 0.5, rho = 0.2, phi.z = 0.99, phi.w = 0.99, 
                               phi.tau = 0.99, lambda = 3, tau.alpha = 1, 
                               tau.beta = 1, nknots = 1, dist = "t")

this.s <- 1
orig.s <- data$y[this.s, 1:(nt - 1)]
lag.s  <- data$y[this.s, 2:nt]
data.s <- cbind(orig.s, lag.s)
plot(data.s)

rank.orig.s <- rank(orig.s) / (nt + 1)
rank.lag.s  <- rank(lag.s) / (nt + 1)
rank.data.s <- cbind(rank.orig.s, rank.lag.s)
plot(rank.data.s, xlim = c(0.90, 1), ylim = c(0.90, 1))

plot(-1 / log(rank.data.s))
chiplot(rank.data.s, which = 2)


# put together matrix to see if partitions are changing
g <- matrix(NA, nt, ns)
for (t in 1:nt) {
  d <- rdist(S, data$knots[, , t])
  for (i in 1:ns) {
    g[t, i] <- which(d[i, ] == min(d[i, ]))
  }
}


cor(data.s3)
chiplot(data.s3)
