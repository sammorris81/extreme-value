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

probs <- c(0.90, 0.95, 0.99)
quantiles.90 <- quantiles.95 <- quantiles.99 <- matrix(NA, dim=c(439, 5))
p.exceed.75 <- array(NA, dim=c(439, 92, 5))

load('./OzoneFull1.RData')
rho <- mean(fit.1$rho)

yp <- fit.1$yp
quantiles.90[, 1] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 1] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 1] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 1] <- mean(yp[, i, t] > 75)
  }
}


load('./OzoneFull2.RData')
yp <- fit.2$yp
quantiles.90[, 2] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 2] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 2] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 2] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull3.RData')
yp <- fit.3$yp
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
quantiles.90[, 5] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 5] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 5] <- apply(yp, 2, quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 5] <- mean(yp[, i, t] > 75)
  }
}

save(quantiles.90, quantiles.95, quantiles.99, p.exceed.75, file="predictions.RData")

# # load("predictions.RData")
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


# z.min <- min(c(quantiles.90[, 4, 1], y[, 4]), na.rm=T)
# z.max <- max(c(quantiles.90[, 4, 1], y[, 4]), na.rm=T)
# quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.90[, 4, 1], nx=25, ny=25, main="Gaussian", col=coltab, zlim=c(z.min, z.max))
# quilt.plot(x=s[(!is.na(y[, 4])), 1], y=s[(!is.na(y[, 4])), 2], z=y[(!is.na(y[, 4])), 4], add=T, add.legend=F, zlim=c(z.min, z.max), col=coltab)
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