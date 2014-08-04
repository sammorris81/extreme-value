rm(list=ls())

probs <- c(0.90, 0.95, 0.99)
quantiles.90 <- quantiles.95 <- quantiles.99 <- array(NA, dim=c(439, 92, 5))
p.exceed.75 <- array(NA, dim=c(439, 92, 5))

load('./OzoneFull1.RData')
yp <- fit.1$yp
quantiles.90[, , 1] <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95[, , 1] <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99[, , 1] <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 1] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull2.RData')
yp <- fit.2$yp
quantiles.90[, , 2] <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95[, , 2] <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99[, , 2] <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 2] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull3.RData')
yp <- fit.3$yp
quantiles.90[, , 3] <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95[, , 3] <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99[, , 3] <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 3] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull4.RData')
yp <- fit.4$yp
quantiles.90[, , 4] <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95[, , 4] <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99[, , 4] <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 4] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull5.RData')
yp <- fit.5$yp
quantiles.90[, , 5] <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95[, , 5] <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99[, , 5] <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 5] <- mean(yp[, i, t] > 75)
  }
}

save(quantiles.90, quantiles.95, quantiles.99, p.exceed.75, file="predictions.RData")

#### Create places for prediciton maps
library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load("predictions.RData")
load('cv-setup-se.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

s1.preds <- seq(1050, 1800, length=25)
s2.preds <- seq(-860, -250, length=25)
s.preds <- expand.grid(s1.preds, s2.preds)
s.preds <- s.preds[(s.preds[, 2] >= (1.33 * s.preds[, 1] - 2815)), ]  # atlantic
s.preds <- s.preds[(s.preds[, 2] < (0.75 * s.preds[, 1] - 1285)), ]  # tennessee
s.preds <- s.preds[(s.preds[, 2] >= (-6.2 * s.preds[, 1] + 5960)), ]  # tennessee

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.90[])
