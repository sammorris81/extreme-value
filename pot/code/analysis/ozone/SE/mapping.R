rm(list=ls())

probs <- c(0.90, 0.95, 0.99)
quantiles <- array(NA, dim=c(5, 439, 3))
p.exceed.75 <- matrix(NA, nrow=5, ncol=439)

load('./OzoneFull1.RData')
yp <- fit.1$yp
for (i in 1:3) {
  quantiles[1, , i] <- apply(yp, 2, quantile, probs=probs[i])
}
# for each site, how much of the posterior predictive is above 75ppb
for (i in 1:439) {
  p.exceed.75[1, i] <- mean(yp[, i, ] > 75)
}ri
load('./OzoneFull2.RData')
yp <- fit.2$yp
for (i in 1:3) {
  quantiles[2, , i] <- apply(yp, 2, quantile, probs=probs[i])
}
# for each site, how much of the posterior predictive is above 75ppb
for (i in 1:439) {
  p.exceed.75[2, i] <- mean(yp[, i, ] > 75)
}

load('./OzoneFull3.RData')
yp <- fit.3$yp
for (i in 1:3) {
  quantiles[3, , i] <- apply(yp, 2, quantile, probs=probs[i])
}
# for each site, how much of the posterior predictive is above 75ppb
for (i in 1:439) {
  p.exceed.75[3, i] <- mean(yp[, i, ] > 75)
}

load('./OzoneFull4.RData')
yp <- fit.4$yp
for (i in 1:3) {
  quantiles[4, , i] <- apply(yp, 2, quantile, probs=probs[i])
}
# for each site, how much of the posterior predictive is above 75ppb
for (i in 1:439) {
  p.exceed.75[4, i] <- mean(yp[, i, ] > 75)
}

load('./OzoneFull5.RData')
yp <- fit.5$yp
for (i in 1:3) {
  quantiles[5, , i] <- apply(yp, 2, quantile, probs=probs[i])
}
# for each site, how much of the posterior predictive is above 75ppb
for (i in 1:439) {
  p.exceed.75[5, i] <- mean(yp[, i, ] > 75)
}

save(quantiles, p.exceed.75, file="predictions.RData")