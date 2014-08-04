rm(list=ls())

probs <- c(0.90, 0.95, 0.99)
quantiles.90 <- quantiles.95 <- quantiles.99 <- array(NA, dim=c(439, 92, 5))
p.exceed.75 <- array(NA, dim=c(439, 92, 5))

load('./OzoneFull1.RData')
yp <- fit.1$yp
quantiles.90 <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95 <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99 <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 1] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull2.RData')
yp <- fit.2$yp
quantiles.90 <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95 <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99 <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 2] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull3.RData')
yp <- fit.3$yp
quantiles.90 <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95 <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99 <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 3] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull4.RData')
yp <- fit.4$yp
quantiles.90 <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95 <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99 <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 4] <- mean(yp[, i, t] > 75)
  }
}

load('./OzoneFull5.RData')
yp <- fit.5$yp
quantiles.90 <- apply(yp, c(2, 3), quantile, probs=0.90)
quantiles.95 <- apply(yp, c(2, 3), quantile, probs=0.95)
quantiles.99 <- apply(yp, c(2, 3), quantile, probs=0.99)
for (i in 1:439) {
  for (t in 1:92) {
    p.exceed.75[i, t, 5] <- mean(yp[, i, t] > 75)
  }
}

save(quantiles, p.exceed.75, file="predictions.RData")

