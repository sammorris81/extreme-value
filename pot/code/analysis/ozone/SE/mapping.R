library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup-se.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

quantiles.90 <- quantiles.95 <- quantiles.99 <- quantiles.999 <- matrix(NA, nrow=439, ncol=6)
p.exceed.75 <- array(NA, dim=c(439, 92, 6))
post.med <- array(NA, dim=c(439, 92, 6))

# Gaussian
load('./OzoneFull1.RData')
yp <- fit.1$yp
post.med[, , 1] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 1] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 1] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 1] <- apply(yp, 2, quantile, probs=0.99)
quantiles.999[, 1] <- apply(yp, 2, quantile, probs=0.999)
p.exceed.75[, , 1] <- apply((yp >= 75), c(2, 3), mean)

# t-1 T-0.9
load('./OzoneFull2.RData')
yp <- fit.2$yp
post.med[, , 2] <- apply(yp, 2, quantile, probs=0.50)
quantiles.90[, 2] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 2] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 2] <- apply(yp, 2, quantile, probs=0.99)
quantiles.999[, 2] <- apply(yp, 2, quantile, probs=0.999)
p.exceed.75[, , 2] <- apply((yp >= 75), c(2, 3), mean)

# skew t-1 T=0.9
load('./OzoneFull3.RData')
yp <- fit.3$yp
post.med[, , 3] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 3] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 3] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 3] <- apply(yp, 2, quantile, probs=0.99)
quantiles.999[, 3] <- apply(yp, 2, quantile, probs=0.999)
p.exceed.75[, , 3] <- apply((yp >= 75), c(2, 3), mean)

# t-3 T=0.9
load('./OzoneFull4.RData')
yp <- fit.4$yp
post.med[, , 4] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 4] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 4] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 4] <- apply(yp, 2, quantile, probs=0.99)
quantiles.999[, 4] <- apply(yp, 2, quantile, probs=0.999)
p.exceed.75[, , 4] <- apply((yp >= 75), c(2, 3), mean)

# skew t-3 T=0.9
load('./OzoneFull5.RData')
yp <- fit.5$yp
post.med[, , 5] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 5] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 5] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 5] <- apply(yp, 2, quantile, probs=0.99)
quantiles.999[, 5] <- apply(yp, 2, quantile, probs=0.999)
p.exceed.75[, , 5] <- apply((yp >= 75), c(2, 3), mean)

# t-1 T=0
load('./OzoneFull6.RData')
yp <- fit.6$yp
post.med[, , 6] <- apply(yp, c(2, 3), quantile, probs=0.50)
quantiles.90[, 6] <- apply(yp, 2, quantile, probs=0.90)
quantiles.95[, 6] <- apply(yp, 2, quantile, probs=0.95)
quantiles.99[, 6] <- apply(yp, 2, quantile, probs=0.99)
quantiles.999[, 6] <- apply(yp, 2, quantile, probs=0.999)
p.exceed.75[, , 6] <- apply((yp >= 75), c(2, 3), mean)

# probability of no exceedance
exceedance.0 <- exceedance.1 <- exceedance.2 <- exceedance.3 <- matrix(0, 439, 6)

for (i in 1:439) {
  exceedance.0[i, 1] <- prod(1 - p.exceed.75[i, , 1])
  exceedance.0[i, 2] <- prod(1 - p.exceed.75[i, , 2])
  exceedance.0[i, 3] <- prod(1 - p.exceed.75[i, , 3])
  exceedance.0[i, 4] <- prod(1 - p.exceed.75[i, , 4])
  exceedance.0[i, 5] <- prod(1 - p.exceed.75[i, , 5])
  exceedance.0[i, 6] <- prod(1 - p.exceed.75[i, , 6])
}

# probability of exactly one exceedance
for(i in 1:439) {
  for (t in 1:92) {
    exceedance.1[i, 1] <- exceedance.1[i, 1] + p.exceed.75[i, t, 1] * prod(1 - p.exceed.75[i, -t, 1])
    exceedance.1[i, 2] <- exceedance.1[i, 2] + p.exceed.75[i, t, 2] * prod(1 - p.exceed.75[i, -t, 2])
    exceedance.1[i, 3] <- exceedance.1[i, 3] + p.exceed.75[i, t, 3] * prod(1 - p.exceed.75[i, -t, 3])
    exceedance.1[i, 4] <- exceedance.1[i, 4] + p.exceed.75[i, t, 4] * prod(1 - p.exceed.75[i, -t, 4])
    exceedance.1[i, 5] <- exceedance.1[i, 5] + p.exceed.75[i, t, 5] * prod(1 - p.exceed.75[i, -t, 5])
    exceedance.1[i, 6] <- exceedance.1[i, 6] + p.exceed.75[i, t, 6] * prod(1 - p.exceed.75[i, -t, 6])
  }
}

# probability of exactly two exceedances
for(i in 1:439) {
  for (t1 in 1:91) { for (t2 in (t1+1):92) {
    exceedance.2[i, 1] <- exceedance.2[i, 1] + prod(p.exceed.75[i, c(t1, t2), 1]) * prod(1 - p.exceed.75[i, -c(t1, t2), 1])
    exceedance.2[i, 2] <- exceedance.2[i, 2] + prod(p.exceed.75[i, c(t1, t2), 2]) * prod(1 - p.exceed.75[i, -c(t1, t2), 2])
    exceedance.2[i, 3] <- exceedance.2[i, 3] + prod(p.exceed.75[i, c(t1, t2), 3]) * prod(1 - p.exceed.75[i, -c(t1, t2), 3])
    exceedance.2[i, 4] <- exceedance.2[i, 4] + prod(p.exceed.75[i, c(t1, t2), 4]) * prod(1 - p.exceed.75[i, -c(t1, t2), 4])
    exceedance.2[i, 5] <- exceedance.2[i, 5] + prod(p.exceed.75[i, c(t1, t2), 5]) * prod(1 - p.exceed.75[i, -c(t1, t2), 5])
    exceedance.2[i, 6] <- exceedance.2[i, 6] + prod(p.exceed.75[i, c(t1, t2), 6]) * prod(1 - p.exceed.75[i, -c(t1, t2), 6])
  } }
}

# # probability of exactly 3 exceedances
# # very very slow. only run on desktop/server
for(i in 1:439) {
  for (t1 in 1:90) { for (t2 in (t1+1):91) { for (t3 in (t2+1):92) {
    # exceedance.3[i, 1] <- exceedance.3[i, 1] + prod(p.exceed.75[i, c(t1, t2, t3), 1]) * prod(1 - p.exceed.75[i, -c(t1, t2, t3), 1])
    # exceedance.3[i, 2] <- exceedance.3[i, 2] + prod(p.exceed.75[i, c(t1, t2, t3), 2]) * prod(1 - p.exceed.75[i, -c(t1, t2, t3), 2])
    # exceedance.3[i, 3] <- exceedance.3[i, 3] + prod(p.exceed.75[i, c(t1, t2, t3), 3]) * prod(1 - p.exceed.75[i, -c(t1, t2, t3), 3])
    # exceedance.3[i, 4] <- exceedance.3[i, 4] + prod(p.exceed.75[i, c(t1, t2, t3), 4]) * prod(1 - p.exceed.75[i, -c(t1, t2, t3), 4])
    # exceedance.3[i, 5] <- exceedance.3[i, 5] + prod(p.exceed.75[i, c(t1, t2, t3), 5]) * prod(1 - p.exceed.75[i, -c(t1, t2, t3), 5])
    # exceedance.3[i, 6] <- exceedance.3[i, 6] + prod(p.exceed.75[i, c(t1, t2, t3), 6]) * prod(1 - p.exceed.75[i, -c(t1, t2, t3), 6])
  } } }
  print(i)
}

exceedance.4 <- 1 - (exceedance.0 + exceedance.1 + exceedance.2 + exceedance.3)

save(post.med, quantiles.90, quantiles.95, quantiles.99, p.exceed.75, exceedance.0, exceedance.1, exceedance.2, exceedance.3, exceedance.4, file="predictions.RData")

# plots of the results
library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('cv-setup-se.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')
load("predictions.RData")
#### Create places for prediciton maps

s1.preds <- seq(1050, 1800, length=25)
s2.preds <- seq(-860, -250, length=25)
s.preds <- expand.grid(s1.preds, s2.preds)
s.preds <- s.preds[(s.preds[, 2] >= (1.33 * s.preds[, 1] - 2815)), ]  # atlantic
s.preds <- s.preds[(s.preds[, 2] < (0.75 * s.preds[, 1] - 1285)), ]  # tennessee
s.preds <- s.preds[(s.preds[, 2] >= (-6.2 * s.preds[, 1] + 5960)), ]  # tennessee

# quantile-95
par(mfrow=c(2, 2))
z.range <- range(quantiles.95[, c(1, 2, 3, 6)])
quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.95[, 1], nx=25, ny=25, zlim=z.range, main="Gaussian")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.95[, 2], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.9")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.95[, 3], nx=25, ny=25, zlim=z.range, main="skew-t, K=1, T=0.9")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.95[, 6], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.0")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.95[, 7], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.9 Site-specific")
lines(l)

# quantile-99
par(mfrow=c(2, 2))
z.range <- range(quantiles.99[, c(1, 2, 3, 6)])
quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.99[, 1], nx=25, ny=25, zlim=z.range, main="Gaussian")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.99[, 2], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.9")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.99[, 3], nx=25, ny=25, zlim=z.range, main="skew-t, K=1, T=0.9")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.99[, 6], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.0")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=quantiles.99[, 7], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.9 Site-specific")
lines(l)

# probability of exceeding standard
par(mfrow=c(2, 2))
z.range <- range(exceedance.4[, c(1, 2, 3, 6)])
quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.4[, 1], nx=25, ny=25, zlim=z.range, main="Gaussian")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.4[, 2], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.9")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.4[, 3], nx=25, ny=25, zlim=z.range, main="skew-t, K=1, T=0.9")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.4[, 6], nx=25, ny=25, zlim=z.range, main="t, K=1, T=0.0")
lines(l)

quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=exceedance.4[, 7], nx=25, ny=25, zlim=z.range, main="skew-t, K=1, T=0.9 Site-specific")
lines(l)