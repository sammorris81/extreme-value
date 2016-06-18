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
phis <- seq(0, 0.98, length = 50)
ext.coef.1 <- ext.coef.3 <- ext.coef.5 <- ext.coef.10 <- rep(0, length(phis))

for (i in 1:length(phis)) {
  # single knot
  set.seed(i * 100)
  phi <- phis[i] 
  ns <- 1
  nt <- 100000
  X <- array(1, dim = c(ns, nt, 1))
  S <- matrix(runif(ns * 2), ns, 2)
  data <- rpotspatTS_cont_lambda(nt = nt, x = X, s = S, beta = 0, gamma = 0.90, 
                                 nu = 0.5, rho = 0.2, phi.z = phi, phi.w = phi, 
                                 phi.tau = phi, lambda = 3, tau.alpha = 1, 
                                 tau.beta = 1, nknots = 1, dist = "t")
  
  cur <- data$y[1:(nt - 1)]
  lag.1 <- data$y[2:nt]
  
  # get f-madogram
  fhat.cur    <- rank(cur) / (length(cur) + 1)
  fhat.lag    <- rank(lag.1) / (length(lag.1) + 1)
  fmado       <- 0.5 * mean(abs(fhat.cur - fhat.lag))
  ext.coef.1[i] <- (1 + 2 * fmado) / (1 - 2 * fmado)
  
  cur <- data$y[1:(nt - 2)]
  lag.3 <- data$y[3:nt]
  
  # get f-madogram
  fhat.cur    <- rank(cur) / (length(cur) + 1)
  fhat.lag    <- rank(lag.3) / (length(lag.3) + 1)
  fmado       <- 0.5 * mean(abs(fhat.cur - fhat.lag))
  ext.coef.3[i] <- (1 + 2 * fmado) / (1 - 2 * fmado)
  
  cur <- data$y[1:(nt - 4)]
  lag.5 <- data$y[5:nt]
  
  # get f-madogram
  fhat.cur    <- rank(cur) / (length(cur) + 1)
  fhat.lag    <- rank(lag.5) / (length(lag.5) + 1)
  fmado       <- 0.5 * mean(abs(fhat.cur - fhat.lag))
  ext.coef.5[i] <- (1 + 2 * fmado) / (1 - 2 * fmado)
  
  cur <- data$y[1:(nt - 9)]
  lag.10 <- data$y[10:nt]
  
  # get f-madogram
  fhat.cur    <- rank(cur) / (length(cur) + 1)
  fhat.lag    <- rank(lag.10) / (length(lag.10) + 1)
  fmado       <- 0.5 * mean(abs(fhat.cur - fhat.lag))
  ext.coef.10[i] <- (1 + 2 * fmado) / (1 - 2 * fmado)
  
  print(paste("i = ", i, sep = ""))
}
save.image(file = "chiphi.RData")

cols <- c("firebrick3", "dodgerblue3", "orange3", "darkolivegreen3")
lambda.1 <- 2 - ext.coef.1
lambda.3 <- 2 - ext.coef.3
lambda.5 <- 2 - ext.coef.5
lambda.10 <- 2 - ext.coef.10
# quartz(width = 4, height = 4)
plot(phis, lambda.1, type = "l",
     xlab = bquote(phi), ylab = bquote(chi), col = cols[1],
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex = 1.5)
lines(phis, lambda.3, col = cols[2])
lines(phis, lambda.5, col = cols[3])
lines(phis, lambda.10, col = cols[4])
legend("topleft", lty=c(1, 1, 1, 1),
       legend=c("Lag 1", "Lag 3", "Lag 5", "Lag 10"),
       col=c("firebrick3", "dodgerblue3", "orange3", "darkolivegreen3"), cex=1.5)
dev.print(device = pdf, file = "plots/chi-phi.pdf")
dev.off()





chiplot(data = data.s1, which = 1, main1 = "site 1", 
        qlim = c(0.005, 0.995), nq = 200, ylim1 = c(0, 1))
chiplot(data = data.s2, which = 1, main1 = "site 2", 
        qlim = c(0.001, 0.999), nq = 200, ylim1 = c(0, 1))
chiplot(data = data.s3, which = 1, main1 = "site 3", qlim = c(0.001, 0.999))
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
