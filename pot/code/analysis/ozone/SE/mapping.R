library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
# load('cv-setup-se.RData')
# source('../../../R/mcmc.R')
# source('../../../R/auxfunctions.R')

# s1.preds <- seq(1050, 1800, length=25)
# s2.preds <- seq(-860, -250, length=25)
# s.preds <- expand.grid(s1.preds, s2.preds)
# s.preds <- s.preds[(s.preds[, 2] >= (1.33 * s.preds[, 1] - 2815)), ]  # atlantic
# s.preds <- s.preds[(s.preds[, 2] < (0.75 * s.preds[, 1] - 1285)), ]  # tennessee
# s.preds <- s.preds[(s.preds[, 2] >= (-6.2 * s.preds[, 1] + 5960)), ]  # tennessee

# s.scale.preds <- matrix(NA, nrow=nrow(s.preds), ncol=ncol(s.preds))
# s.scale.preds[,1] <- (s.preds[,1] - range(s[,1])[1])/(range(s[,1])[2] - range(s[,1])[1])
# s.scale.preds[,2] <- (s.preds[,2] - range(s[,2])[1])/(range(s[,2])[2] - range(s[,2])[1])

# X <- X[, , c(1, 2, 3)]
# X.preds <- array(1, dim=c(nrow(s.preds), nt, 3))
# for (t in 1:nt) {
  # X.preds[, t, 2] <- s.scale.preds[, 1]
  # X.preds[, t, 3] <- s.scale.preds[, 2]
# }

probs <- c(0.50, 0.90, 0.95, 0.99)
quantiles.90 <- quantiles.95 <- quantiles.99 <- matrix(NA, nrow=439, ncol=5)
p.exceed.75 <- array(NA, dim=c(439, 92, 5))
post.med <- array(NA, dim=c(439, 92, 5))

tau <- matrix(NA, 92, 3)
beta <- matrix(NA, 3, 3)
rho <- rep(NA, 3)
nu <- rep(NA, 3)
alpha <- rep(NA, 3)


load('./OzoneFull1.RData')
tau[, 1] <- apply(fit.1$tau, 2, mean)
beta[, 1] <- apply(fit.1$beta, 2, mean)
rho[1] <- mean(fit.1$rho)
nu[1] <- mean(fit.1$nu)
alpha[1] <- mean(fit.1$alpha)
yp <- fit.1$yp
post.med[, , 1] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 1] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 1] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 1] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 1] <- mean(yp[, i, t] > 75)
  }
}


load('./OzoneFull2.RData')
tau[, 2] <- apply(fit.2$tau, 2, mean)
beta[, 2] <- apply(fit.2$beta, 2, mean)
rho[2] <- mean(fit.2$rho)
nu[2] <- mean(fit.2$nu)
alpha[2] <- mean(fit.2$alpha)
yp <- fit.2$yp
post.med[, , 2] <- apply(yp, 2, quantile, probs=0.50)
quantiles.90[, 2] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 2] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 2] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 2] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull3.RData')
tau[, 3] <- apply(fit.3$tau, 2, mean)
beta[, 3] <- apply(fit.3$beta, 2, mean)
rho[3] <- mean(fit.3$rho)
nu[3] <- mean(fit.3$nu)
alpha[3] <- mean(fit.3$alpha)
z.alpha <- mean(fit.3$z.alpha)
z <- apply(fit.3$z, c(2, 3), mean)
yp <- fit.3$yp
post.med[, , 3] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 3] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 3] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 3] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 3] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull4.RData')
yp <- fit.4$yp
post.med[, , 4] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 4] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 4] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 4] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 4] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull5.RData')
yp <- fit.5$yp
post.med[, , 5] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 5] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 5] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 5] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 5] <- mean(yp[, i, t] > 75)
  }
}

save(post.med, quantiles.90, quantiles.95, quantiles.99, p.exceed.75, tau, beta, rho, nu, alpha, z.alpha, file="predictions.RData")

# # plot monitoring station ozone locations
# plot(s, type="p", xlab="", ylab="")
# lines(l)

# # plot ozone for days 5 and 34
# par(mfrow=c(1, 2))
# zlim=range(y[, c(5, 34)], na.rm=T)
# quilt.plot(x=s[, 1], y=s[, 2], z=y[, 5], nx=40, ny=40, zlim=zlim, main="Day 5", xlab="", ylab="")
# lines(l)
# quilt.plot(x=s[, 1], y=s[, 2], z=y[, 34], nx=40, ny=40, zlim=zlim, main="Day 34", xlab="", ylab="")
# lines(l)

# load("predictions.RData")
# #### Create places for prediciton maps
# s1.preds <- seq(1050, 1800, length=25)
# s2.preds <- seq(-860, -250, length=25)
# s.preds <- expand.grid(s1.preds, s2.preds)
# s.preds <- s.preds[(s.preds[, 2] >= (1.33 * s.preds[, 1] - 2815)), ]  # atlantic
# s.preds <- s.preds[(s.preds[, 2] < (0.75 * s.preds[, 1] - 1285)), ]  # tennessee
# s.preds <- s.preds[(s.preds[, 2] >= (-6.2 * s.preds[, 1] + 5960)), ]  # tennessee

# s.scale.preds <- matrix(NA, nrow=nrow(s.preds), ncol=ncol(s.preds))
# s.scale.preds[,1] <- (s.preds[,1] - range(s[,1])[1])/(range(s[,1])[2] - range(s[,1])[1])
# s.scale.preds[,2] <- (s.preds[,2] - range(s[,2])[1])/(range(s[,2])[2] - range(s[,2])[1])

# X <- X[, , c(1, 2, 3)]
# X.preds <- array(1, dim=c(nrow(s.preds), nt, 3))
# for (t in 1:nt) {
  # X.preds[, t, 2] <- s.scale.preds[, 1]
  # X.preds[, t, 3] <- s.scale.preds[, 2]
# }


# exceedance.1 <- exceedance.2 <- exceedance.3 <- exceedance.4 <- exceedance.5 <- rep(NA, 439)
# for (i in 1:439) {
  # exceedance.1[i] <- 1 - prod(1 - p.exceed.75[i, , 1])
  # exceedance.2[i] <- 1 - prod(1 - p.exceed.75[i, , 2])
  # exceedance.3[i] <- 1 - prod(1 - p.exceed.75[i, , 3])
  # exceedance.4[i] <- 1 - prod(1 - p.exceed.75[i, , 4])
  # exceedance.5[i] <- 1 - prod(1 - p.exceed.75[i, , 5])
# }

# samp.p.exceed.75 <- rep(NA, 85)
# for (i in 1:85) {
  # samp.p.exceed.75[i] <- mean(y[i, ] > 75, na.rm=T)
# }

# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=post.med[, 5, 1], nx=25, ny=25)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=post.med[, 34, 1], nx=25, ny=25)

# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=post.med[, 5, 4], nx=25, ny=25)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=post.med[, 34, 4], nx=25, ny=25)

# coltab <- tim.colors(256)
# z.min <- min(samp.p.exceed.75, exceedance.1)
# z.max <- max(samp.p.exceed.75, exceedance.1)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.1, nx=25, ny=25, main="Gaussian", col=coltab, zlim=c(z.min, z.max))
# quilt.plot(x=s[, 1], y=s[, 2], z=samp.p.exceed.75, main="sample probability", add=T, add.legend=F, zlim=c(z.min, z.max), col=coltab)
# lines(l)

# z.min <- min(samp.p.exceed.75, exceedance.2)
# z.max <- max(samp.p.exceed.75, exceedance.2)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.2, nx=25, ny=25, main="t, K=1, T=0.90", col=coltab, zlim=c(z.min, z.max))
# quilt.plot(x=s[, 1], y=s[, 2], z=samp.p.exceed.75, main="sample probability", add=T, add.legend=F, zlim=c(z.min, z.max), col=coltab)
# lines(l)

# z.min <- min(samp.p.exceed.75, exceedance.3)
# z.max <- max(samp.p.exceed.75, exceedance.3)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.3, nx=25, ny=25, main="skew-t, K=1, T=0.90", col=coltab, zlim=c(z.min, z.max))
# quilt.plot(x=s[, 1], y=s[, 2], z=samp.p.exceed.75, main="sample probability", add=T, add.legend=F, zlim=c(z.min, z.max), col=coltab)
# lines(l)

# z.min <- min(samp.p.exceed.75, exceedance.4)
# z.max <- max(samp.p.exceed.75, exceedance.4)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.4, nx=25, ny=25, main="t, K=3, T=0.90", col=coltab, zlim=c(z.min, z.max))
# quilt.plot(x=s[, 1], y=s[, 2], z=samp.p.exceed.75, main="sample probability", add=T, add.legend=F, zlim=c(z.min, z.max), col=coltab)
# lines(l)


# z.min <- min(samp.p.exceed.75, exceedance.5)
# z.max <- max(samp.p.exceed.75, exceedance.5)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.5, nx=25, ny=25, main="skew-t, K=3, T=0.90", col=coltab, zlim=c(z.min, z.max))
# quilt.plot(x=s[, 1], y=s[, 2], z=samp.p.exceed.75, main="sample probability", add=T, add.legend=F, zlim=c(z.min, z.max), col=coltab)
# lines(l)


# z.min <- min(quantiles.90[, 1], na.rm=T)
# z.max <- max(quantiles.90[, 1], na.rm=T)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.90[, 1], nx=25, ny=25, main="Gaussian", col=coltab, zlim=c(z.min, z.max))
# lines(l)

# z.min <- min(quantiles.95[, 2], na.rm=T)
# z.max <- max(quantiles.95[, 2], na.rm=T)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.95[, 2], nx=25, ny=25, main="t, K=1, T=0.90", col=coltab, zlim=c(z.min, z.max))
# lines(l)


# s1.preds <- seq(1050, 1800, length=30)
# s2.preds <- seq(-860, -250, length=30)
# s.preds <- expand.grid(s1.preds, s2.preds)
# knots <- cbind(runif(5, 1050, 1800), runif(5, -860, -250))
# g <- mem(s.preds, knots)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=g, nx=30, ny=30, add.legend=F)
# text(knots[1, 1], knots[1, 2], "1", cex=4)
# text(knots[2, 1], knots[2, 2], "2", cex=4)
# text(knots[3, 1], knots[3, 2], "3", cex=4)
# text(knots[4, 1], knots[4, 2], "4", cex=4)
# text(knots[5, 1], knots[5, 2], "5", cex=4)
# lines(l)

# what if I took the posterior mean for each of the model parameters (gaussian and t-1 T=0.9 and skew-t-1 T=0.90 and just generated the data)