load("us-all-setup.RData")
source("../../../R/auxfunctions.R")
settings <- read.csv("settings.csv")

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
thresholds <- quantile(Y, probs=probs, na.rm=T)
nsets <- 2 # Number of cv sets
nbetas <- 2 # number of betas

quant.score <- array(NA, dim=c(length(probs), nsets, 50))
brier.score <- array(NA, dim=c(length(thresholds), nsets, 50))

beta.0 <- array(NA, dim=c(5000, nsets, 50))
beta.1 <- array(NA, dim=c(5000, nsets, 50))

phi.z <- array(NA, dim=c(5000, nsets, 24))
phi.w <- array(NA, dim=c(5000, nsets, 24))
phi.tau <- array(NA, dim=c(5000, nsets, 24))

load("us-all-results.RData")
done <- c(43:50)
quant.score <- savelist[[1]]
brier.score <- savelist[[2]]
beta.0 <- savelist[[3]]
beta.1 <- savelist[[4]]
probs <- savelist[[5]]
thresholds <- savelist[[6]]

for (i in 1:50) {
  file <- paste("us-all-", i, ".RData", sep="")
  cat("start file", file, "\n")
  if (i %in% done) {
    load(file)
    for (d in 1:nsets) {
      fit.d <- fit[[d]]
      val.idx <- cv.lst[[d]]
      validate <- Y[val.idx, ]
      pred.d <- fit.d$yp[, , ]
      if (i == 26) {
        validate <- t(validate)
        pred.d <- fit.d$yp[(25001:30000), , ]
      }
      quant.score[, d, i] <- QuantScore(pred.d, probs, validate)
      brier.score[, d, i] <- BrierScore(pred.d, thresholds, validate)
      if (i != 26) {
        beta.0[, d, i] <- fit.d$beta[, 1]
        beta.1[, d, i] <- fit.d$beta[, 2]
      }
    }
  }
  cat("finish file", file, "\n")
}

savelist <- list(quant.score, brier.score,
                 beta.0, beta.1,
                 probs, thresholds)

save(savelist, file="us-all-results.RData")

rm(list=ls())
load("us-all-setup.RData")
source("../../../R/auxfunctions.R")
load("us-all-results.RData")
settings <- read.csv("settings.csv")

quant.score <- savelist[[1]]
brier.score <- savelist[[2]]
beta.0 <- savelist[[3]]
beta.1 <- savelist[[4]]
probs <- savelist[[5]]
thresholds <- savelist[[6]]

quant.score.mean <- matrix(NA, 50, length(probs))
brier.score.mean <- matrix(NA, 50, length(thresholds))

quant.score.se <- matrix(NA, 50, length(probs))
brier.score.se <- matrix(NA, 50, length(thresholds))

done <- c(1:25, 27:50)
for (i in 1:50) {
  if (i %in% done) {
    quant.score.mean[i, ] <- apply(quant.score[, , i], 1, mean)
    quant.score.se[i, ] <- apply(quant.score[, , i], 1, sd) / sqrt(2)
    brier.score.mean[i, ] <- apply(brier.score[, , i], 1, mean)
    brier.score.se[i, ] <- apply(brier.score[, , i], 1, sd) / sqrt(2)
  }
}

quant.score.mean[c(1:10, 13, 14, 17:19, 21, 23, 25, 26), ]
for (i in 1:length(thresholds)) {
 print(which(quant.score.mean[, i] == min(quant.score.mean[, i], na.rm=T)))
}

for (i in 1:length(thresholds)) {
  print(which(brier.score.mean[, i] == min(brier.score.mean[, i], na.rm=T)))
}

bs.mean.ref.gau <- matrix(NA, nrow=49, ncol=11)
qs.mean.ref.gau <- matrix(NA, nrow=49, ncol=11)
for (i in 1:49) {
  bs.mean.ref.gau[i, ] <- brier.score.mean[(i + 1), ] / brier.score.mean[1, ]
  qs.mean.ref.gau[i, ] <- quant.score.mean[(i + 1), ] / quant.score.mean[1, ]
}

bg <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1", "gray80")
col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4", "gray16")

# one knot
quartz(width=15, height=12)
par(mfrow=c(3, 2), mar=c(5.1, 5.1, 4.1, 2.1))
xplot <- probs[6:11]
plot(xplot, bs.mean.ref.gau[1, 6:11], type="b", ylim=c(0.8, 1.2), pch=21,
     col=col[1], bg=bg[1], ylab="brier score", xlab="sample quantiles",
     main="Brier scores (K = 1)", cex.lab=2, cex.axis=2, cex.main=2, cex=1.7)
lines(xplot, bs.mean.ref.gau[2, 6:11], type="b", pch=21, col=col[2], bg=bg[2], cex=1.7)   # T = 50
lines(xplot, bs.mean.ref.gau[3, 6:11], type="b", pch=21, col=col[3], bg=bg[3], cex=1.7)   # T = 75
abline(h=1, lty=2)
lines(xplot, bs.mean.ref.gau[4, 6:11], type="b", pch=21, col=col[4], bg=bg[4], cex=1.7)   # T = 90
lines(xplot, bs.mean.ref.gau[17, 6:11], type="b", pch=24, col=col[3], bg=bg[3], cex=1.7)  # T = 75, no skew
lines(xplot, bs.mean.ref.gau[18, 6:11], type="b", pch=24, col=col[4], bg=bg[4], cex=1.7)  # T = 90, no skew
# plot(bs.mean.ref.gau[26, ], type="b")  # T = 0, time series
# lines(bs.mean.ref.gau[27, ], type="b")  # T = 50, time series
# lines(bs.mean.ref.gau[28, ], type="b")  # T = 75, time series
# lines(bs.mean.ref.gau[29, ], type="b")  # T = 90, time series
# lines(bs.mean.ref.gau[42, ], type="b")  # T = 75, time series, no skew
# lines(bs.mean.ref.gau[43, ], type="b")  # T = 90, time series, no skew

# 5 knots
plot(xplot, bs.mean.ref.gau[5, 6:11], type="b", ylim=c(0.8, 1.2), pch=21,
     col=col[1], bg=bg[1], ylab="brier score", xlab="sample quantiles",
     main="Brier scores (K = 5)", cex.lab=2, cex.axis=2, cex.main=2, cex=1.7)
lines(xplot, bs.mean.ref.gau[6, 6:11], type="b", pch=21, col=col[2], bg=bg[2], cex=1.7)  # unusually high
lines(xplot, bs.mean.ref.gau[7, 6:11], type="b", pch=21, col=col[3], bg=bg[3], cex=1.7)
abline(h=1, lty=2)
lines(xplot, bs.mean.ref.gau[8, 6:11], type="b", pch=21, col=col[4], bg=bg[4], cex=1.7)
lines(xplot, bs.mean.ref.gau[19, 6:11], type="b", pch=24, col=col[3], bg=bg[3], cex=1.7)  # T = 75, no skew
lines(xplot, bs.mean.ref.gau[20, 6:11], type="b", pch=24, col=col[4], bg=bg[4], cex=1.7)  # T = 90, no skew
# lines(bs.mean.ref.gau[30, ], type="b")  # T = 0, time series
# lines(bs.mean.ref.gau[31, ], type="b")  # T = 50, time series
# lines(bs.mean.ref.gau[32, ], type="b")  # T = 75, time series
# lines(bs.mean.ref.gau[33, ], type="b")  # T = 90, time series
# lines(bs.mean.ref.gau[44, ], type="b")  # T = 75, time series, no skew
# lines(bs.mean.ref.gau[45, ], type="b")  # T = 90, time series, no skew

# 10 knots
plot(xplot, bs.mean.ref.gau[9, 6:11], type="b", ylim=c(0.8, 1.2), pch=21,
     col=col[1], bg=bg[1], ylab="brier score", xlab="sample quantiles",
     main="Brier scores (K = 10)", cex.lab=2, cex.axis=2, cex.main=2, cex=1.7)
lines(xplot, bs.mean.ref.gau[10, 6:11], type="b", pch=21, col=col[2], bg=bg[2], cex=1.7)
lines(xplot, bs.mean.ref.gau[11, 6:11], type="b", pch=21, col=col[3], bg=bg[3], cex=1.7)
abline(h=1, lty=2)
lines(xplot, bs.mean.ref.gau[12, 6:11], type="b", pch=21, col=col[4], bg=bg[4], cex=1.7)
lines(xplot, bs.mean.ref.gau[21, 6:11], type="b", pch=24, col=col[3], bg=bg[3], cex=1.7)  # T = 75, no skew
lines(xplot, bs.mean.ref.gau[22, 6:11], type="b", pch=24, col=col[4], bg=bg[4], cex=1.7)  # T = 90, no skew
# lines(bs.mean.ref.gau[34, ], type="b")  # T = 0, time series
# lines(bs.mean.ref.gau[35, ], type="b")  # T = 50, time series
# lines(bs.mean.ref.gau[36, ], type="b")  # T = 75, time series
# lines(bs.mean.ref.gau[37, ], type="b")  # T = 90, time series
# lines(bs.mean.ref.gau[46, ], type="b")  # T = 75, time series, no skew
# lines(bs.mean.ref.gau[47, ], type="b")  # T = 90, time series, no skew

# 15 knots
plot(xplot, bs.mean.ref.gau[13, 6:11], type="b", ylim=c(0.8, 1.2), pch=21,
  col=col[1], bg=bg[1], ylab="brier score", xlab="sample quantiles",
  main="Brier scores (K = 15)", cex.lab=2, cex.axis=2, cex.main=2, cex=1.7)
lines(xplot, bs.mean.ref.gau[14, 6:11], type="b", pch=21, col=col[2], bg=bg[2], cex=1.7)  # unusually high
lines(xplot, bs.mean.ref.gau[15, 6:11], type="b", pch=21, col=col[3], bg=bg[3], cex=1.7)
abline(h=1, lty=2)
lines(xplot, bs.mean.ref.gau[16, 6:11], type="b", pch=21, col=col[4], bg=bg[4], cex=1.7)
lines(xplot, bs.mean.ref.gau[23, 6:11], type="b", pch=24, col=col[3], bg=bg[3], cex=1.7)  # T = 75, no skew
lines(xplot, bs.mean.ref.gau[24, 6:11], type="b", pch=24, col=col[4], bg=bg[4], cex=1.7)  # T = 90, no skew
# lines(bs.mean.ref.gau[38, ], type="b")  # T = 0, time series
# lines(bs.mean.ref.gau[39, ], type="b")  # T = 50, time series
# lines(bs.mean.ref.gau[40, ], type="b")  # T = 75, time series
# lines(bs.mean.ref.gau[41, ], type="b")  # T = 90, time series
# lines(bs.mean.ref.gau[48, ], type="b")  # T = 75, time series, no skew
# lines(bs.mean.ref.gau[49, ], type="b")  # T = 90, time series, no skew

# max-stable
plot(xplot, bs.mean.ref.gau[25, 6:11], type="b", ylim=c(0.8, 1.2), pch=23,
     col=col[5], bg=bg[5], ylab="brier score", xlab="sample quantiles",
     main="Brier scores (Max stable)", cex.lab=2, cex.axis=2, cex.main=2, cex=1.7)
abline(h=1, lty=2)

plot(xplot, type="n", axes=F, ylab="", xlab="")
legend("center", legend=c("Skew-t, T=0", "Skew-t, T=50", "Skew-t, T=75", "Skew-t, T=85", "T, T=75", "T, T=85"),
       col=c(col[1:4], col[3:4]), pch=c(rep(21, 4), rep(24, 2)), pt.bg=c(bg[1:4], bg[3:4]), cex=2.4, box.lty=0)



# best performing partition settings
plot(xplot, bs.mean.ref.gau[1, these.probs], type="b", ylim=c(0.755, 0.843), pch=21,
     col=col[1], bg=bg[1], ylab="Relative Brier score", xlab="Threshold quantile", lty=1,
     # main="Select ozone Brier scores",
     cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex=1.3)
lines(xplot, bs.mean.ref.gau[5, these.probs], type="b", pch=21, col=col[2], bg=bg[2],
      cex=1.3, lty=1)
lines(xplot, bs.mean.ref.gau[9, these.probs], type="b", pch=21, col=col[3], bg=bg[3],
      cex=1.3, lty=1)
lines(xplot, bs.mean.ref.gau[13, these.probs], type="b", pch=21, col=col[4], bg=bg[4],
      cex=1.3, lty=1)
abline(h=1, lty=2)
legend("topleft", legend=c("Skew-t, K=1, T=0", "Skew-t, K=5, T=0", "Skew-t, K=10, T=0",
       "Skew-t, K=15, T=0"),
       pch=c(21, 21, 21, 21), lty=c(1, 1, 1, 1), cex=1.3,
       pt.bg=c(bg[1], bg[2], bg[3], bg[4]),
       col=c(col[1], col[2], col[3], col[4]))

bg <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1", "gray80")
col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4", "gray16")
these.probs <- 1:10
xplot <- probs[these.probs]
plot(xplot, bs.mean.ref.gau[32, these.probs], type="b", ylim=c(0.755, 0.843), pch=21,
     col=col[1], bg=bg[1], ylab="Relative Brier score", xlab="Threshold quantile", lty=1,
     # main="Select ozone Brier scores",
     cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex=1.3)
lines(xplot, bs.mean.ref.gau[33, these.probs], type="b", pch=21, col=col[2], bg=bg[2],
      cex=1.3, lty=1)
lines(xplot, bs.mean.ref.gau[34, these.probs], type="b", pch=21, col=col[3], bg=bg[3],
      cex=1.3, lty=1)
lines(xplot, bs.mean.ref.gau[35, these.probs], type="b", pch=21, col=col[4], bg=bg[4],
      cex=1.3, lty=1)
abline(h=1, lty=2)
legend("topleft", legend=c("Skew-t, K=6, T=0", "Skew-t, K=7, T=0", "Skew-t, K=8, T=0",
       "Skew-t, K=9, T=0"),
       pch=c(21, 21, 21, 21), lty=c(1, 1, 1, 1), cex=1.3,
       pt.bg=c(bg[1], bg[2], bg[3], bg[4]),
       col=c(col[1], col[2], col[3], col[4]))



wilcox.results.gau <- matrix(NA, nrow=length(probs), ncol=4)
part.best <- c(6, 33, 34, 35)
for(i in 1:length(probs)) { for (k in 1:4) {
  wilcox.results.gau[i, k] <- wilcox.test(brier.score[i, , 1], brier.score[i, , 6])$p.value
}}

bg <- c("firebrick1", "dodgerblue1")
col <- c("firebrick4", "dodgerblue4")

# T = 0
xplot <- probs
plot(xplot, bs.mean.ref.gau[1, ], type="b", ylim=c(0.8, 1.05), col=col[1], bg=bg[1], pch=19)
lines(xplot, bs.mean.ref.gau[5, ], type="b", col=col[1], bg=bg[1], pch=19)   # 5 knots
lines(xplot, bs.mean.ref.gau[9, ], type="b", col=col[1], bg=bg[1], pch=19)   # 10 knots
lines(xplot, bs.mean.ref.gau[13, ], type="b", col=col[1], bg=bg[1], pch=19)  # 15 knots
lines(xplot, bs.mean.ref.gau[26, ], type="b", col=col[2], bg=bg[2], pch=19)  # 1 knot, time series
lines(xplot, bs.mean.ref.gau[30, ], type="b", col=col[2], bg=bg[2], pch=19)  # 5 knots, time series
lines(xplot, bs.mean.ref.gau[34, ], type="b", col=col[2], bg=bg[2], pch=19)  # 10 knots, time series
lines(xplot, bs.mean.ref.gau[38, ], type="b", col=col[2], bg=bg[2], pch=19)  # 15 knots, time series

# T = 50
plot(xplot, bs.mean.ref.gau[2, ], type="b", ylim=c(0.75, 1.05), col=col[1], bg=bg[1], pch=19)
lines(xplot, bs.mean.ref.gau[6, ], type="b", col=col[1], bg=bg[1], pch=19)   # 5 knots
lines(xplot, bs.mean.ref.gau[10, ], type="b", col=col[1], bg=bg[1], pch=19)  # 10 knots
lines(xplot, bs.mean.ref.gau[14, ], type="b", col=col[1], bg=bg[1], pch=19)  # 15 knots
lines(xplot, bs.mean.ref.gau[27, ], type="b", col=col[2], bg=bg[2], pch=19)  # 1 knot, time series
lines(xplot, bs.mean.ref.gau[31, ], type="b", col=col[2], bg=bg[2], pch=19)  # 5 knots, time series
lines(xplot, bs.mean.ref.gau[35, ], type="b", col=col[2], bg=bg[2], pch=19)  # 10 knots, time series
lines(xplot, bs.mean.ref.gau[39, ], type="b", col=col[2], bg=bg[2], pch=19)  # 15 knots, time series

# T = 75
xplot <- probs[6:11]
plot(xplot, bs.mean.ref.gau[3, c(6:11)], type="b", ylim=c(0.5, 1.5), col=col[1], bg=bg[1], pch=19)
lines(xplot, bs.mean.ref.gau[7, c(6:11)], type="b", col=col[1], bg=bg[1], pch=19)   # 5 knots
lines(xplot, bs.mean.ref.gau[11, c(6:11)], type="b", col=col[1], bg=bg[1], pch=19)  # 10 knots
lines(xplot, bs.mean.ref.gau[15, c(6:11)], type="b", col=col[1], bg=bg[1], pch=19)  # 15 knots
lines(xplot, bs.mean.ref.gau[17, c(6:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 1 knot, no skew
lines(xplot, bs.mean.ref.gau[19, c(6:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 5 knots, no skew
lines(xplot, bs.mean.ref.gau[21, c(6:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 10 knots, no skew
lines(xplot, bs.mean.ref.gau[23, c(6:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 15 knots, no skew
lines(xplot, bs.mean.ref.gau[28, c(6:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 1 knot, time series
lines(xplot, bs.mean.ref.gau[32, c(6:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 5 knots, time series
lines(xplot, bs.mean.ref.gau[36, c(6:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 10 knots, time series
lines(xplot, bs.mean.ref.gau[40, c(6:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 15 knots, time series
lines(xplot, bs.mean.ref.gau[42, c(6:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 1 knot, no skew, time series
lines(xplot, bs.mean.ref.gau[44, c(6:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 5 knots, no skew, time series
lines(xplot, bs.mean.ref.gau[46, c(6:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 10 knots, no skew, time series
lines(xplot, bs.mean.ref.gau[48, c(6:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 15 knots, no skew, time series

# T = 90
xplot <- probs[10:11]
plot(xplot, bs.mean.ref.gau[4, c(10:11)], type="b", ylim=c(0.5, 4), col=col[1], bg=bg[1], pch=19)
lines(xplot, bs.mean.ref.gau[8, c(10:11)], type="b", col=col[1], bg=bg[1], pch=19)   # 5 knots
lines(xplot, bs.mean.ref.gau[12, c(10:11)], type="b", col=col[1], bg=bg[1], pch=19)  # 10 knots
lines(xplot, bs.mean.ref.gau[16, c(10:11)], type="b", col=col[1], bg=bg[1], pch=19)  # 15 knots
lines(xplot, bs.mean.ref.gau[18, c(10:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 1 knot, no skew
lines(xplot, bs.mean.ref.gau[20, c(10:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 5 knots, no skew
lines(xplot, bs.mean.ref.gau[22, c(10:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 10 knots, no skew
lines(xplot, bs.mean.ref.gau[24, c(10:11)], type="b", col=col[1], bg=bg[1], pch=15)  # 15 knots, no skew
lines(xplot, bs.mean.ref.gau[29, c(10:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 1 knot, time series
lines(xplot, bs.mean.ref.gau[33, c(10:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 5 knots, time series
lines(xplot, bs.mean.ref.gau[37, c(10:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 10 knots, time series
lines(xplot, bs.mean.ref.gau[41, c(10:11)], type="b", col=col[2], bg=bg[2], pch=19)  # 15 knots, time series
lines(xplot, bs.mean.ref.gau[43, c(10:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 1 knot, no skew, time series
lines(xplot, bs.mean.ref.gau[45, c(10:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 5 knots, no skew, time series
lines(xplot, bs.mean.ref.gau[47, c(10:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 10 knots, no skew, time series
lines(xplot, bs.mean.ref.gau[49, c(10:11)], type="b", col=col[2], bg=bg[2], pch=15)  # 15 knots, no skew, time series

library(fields)
xplot <- c(1:4)
yplot <- c(1:4)
z <- matrix(quant.score.mean[c(2:17),6], nrow=4, ncol=4, byrow=F)
zlim=range(z)
image.plot(x=c(1:4), y=c(1:4), z=z, axes=F, xlab="", ylab="", main="Quantile score for q(0.95)", zlim=zlim)
text(c(row(z)), c(col(z)), "")  # not sure why we need this, but without it the labels don't turn on...
axis(side=1, at=c(1:4), labels=c("T=0", "T=50", "T=75", "T=90"), line=0.3, lwd=0)
axis(side=2, at=c(1:4), labels=c("K=1", "K=5", "K=10", "K=15"), line=0.3, las=2, lwd=0)
abline(h=c(1.5, 2.5, 3.5))
abline(v=c(1.5, 2.5, 3.5))

include.qscore <- c(1, 10, 18)
pch <- c(1, 2, 5)
plot(probs, quant.score.mean[1, ], lty=1, type="b", ylim=c(min(quant.score.mean[c(1, 10, 18), ]), max(quant.score.mean[c(1, 10, 18)])), main="Quantile Scores", xlab="quantile", ylab="score")
for (i in 2:3) {
  lines(probs, quant.score.mean[include.qscore[i], ], lty=1, type="b", pch=pch[i])
}
legend("topright", pch=c(21, 24, 23), lty=1, legend=c("Gaussian", "skew-t, K=10, T=q(0)", "Max-stable"), pt.bg="white")

round(quant.score.mean[c(2,6,10,14),c(6, 9:12)], 4)
round(quant.score.se[c(10:13,15:19),c(6, 9:12)], 4)
round(brier.score.mean[c(10:19),c(6, 9:12)]*1000, 3)
round(brier.score.se[c(10:19),c(6, 9:12)]*1000, 3)

round(quant.score.mean, 4)
round(brier.score.mean*1000, 4)

# posterior prediction maps
rm(list=ls())
library(fields)
library(SpatialTools)
load('us-all-setup.RData')
# select every third CMAQ value
S.p <- expand.grid(x, y)
keep.these <- which((S.p[, 1] > -2.3) & (S.p[, 1] < 2.4) & S.p[, 2] > -1.6 & S.p[, 2] < 1.3)
CMAQ.p <- CMAQ[keep.these, ]
S.p    <- S.p[keep.these, ]
nx <- length(unique(S.p[, 1]))
ny <- length(unique(S.p[, 2]))

#### Thin the rows and columns by 3
# First, figure out what the x and y values are for the rows
# and columns we should keep
unique.x <- unique(S.p[, 1])
keep.x <- unique.x[seq(1, length(unique.x), by=10)]
nx <- length(keep.x)
unique.y <- unique(S.p[, 2])
keep.y <- unique.y[seq(1, length(unique(S.p[, 2])), by=10)]
ny <- length(keep.y)
keep.these <- which((S.p[, 1] %in% keep.x) & (S.p[, 2] %in% keep.y))

# Now select the subset from the original covariate and location information
CMAQ.p <- CMAQ.p[keep.these, ]
S.p    <- S.p[keep.these, ]

# load results
threshold <- 75
load('us-all-full-1.RData')
np <- dim(fit$yp)[2]
gaus.95 <- apply(fit$yp, c(2, 3), quantile, probs=0.95)
gaus.p.below <- matrix(0, np, nt)
for (i in 1:np) {
  gaus.p.exceed <-
}
load('us-all-full-2.RData')
t1.yp <- fit$yp
load('us-all-full-33.RData')
t6.yp <- fit$yp
rm(fit)

# zlim=c(0, 122)
# quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(CMAQ.p[, 5]), zlim=zlim, nx=nx, ny=ny)
# lines(borders/1000)

# posterior maps: 95th quantile
gaus.95 <- apply(gaus.yp, c(2, 3), quantile, probs=0.95)

# probability of exceeding 75ppb at least once



load("us-all-setup.RData")
source("../../../R/auxfunctions.R")

load("us-all-10.RData")

par(mfrow=c(5, 5))
days1 <- c(1, 4, 7, 10, 13)
days2 <- c(16, 19, 22, 25, 28)

# plot tau
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

# plot z
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }


load("us-all-6.RData")
days1 <- c(1, 4, 7, 10, 13)
days2 <- c(16, 19, 22, 25, 28)
par(mfrow=c(5, 5))
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

par(mfrow=c(2, 5))
plot(fit[[1]]$alpha, type="l")
plot(fit[[1]]$tau.alpha, type="l")
plot(fit[[1]]$tau.beta, type="l")
plot(fit[[1]]$rho, type="l")
plot(fit[[1]]$nu, type="l")
plot(fit[[1]]$beta[, 1], type="l")
plot(fit[[1]]$beta[, 2], type="l")

par(mfrow=c(5, 5))
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[2]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[2]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[2]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[2]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

par(mfrow=c(2, 5))
plot(fit[[2]]$alpha, type="l")
plot(fit[[2]]$tau.alpha, type="l")
plot(fit[[2]]$tau.beta, type="l")
plot(fit[[2]]$rho, type="l")
plot(fit[[2]]$nu, type="l")
plot(fit[[2]]$beta[, 1], type="l")
plot(fit[[2]]$beta[, 2], type="l")

load("us-all-31.RData")
days1 <- c(1, 4, 7, 10, 13)
days2 <- c(16, 19, 22, 25, 28)
par(mfrow=c(5, 5))
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

par(mfrow=c(2, 5))
plot(fit[[1]]$alpha, type="l")
plot(fit[[1]]$tau.alpha, type="l")
plot(fit[[1]]$tau.beta, type="l")
plot(fit[[1]]$rho, type="l")
plot(fit[[1]]$nu, type="l")
plot(fit[[1]]$phi.z, type="l")
plot(fit[[1]]$phi.tau, type="l")
plot(fit[[1]]$phi.w, type="l")
plot(fit[[1]]$beta[, 1], type="l")
plot(fit[[1]]$beta[, 2], type="l")

par(mfrow=c(5, 5))
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[2]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[2]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[2]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[2]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

par(mfrow=c(2, 5))
plot(fit[[2]]$alpha, type="l")
plot(fit[[2]]$tau.alpha, type="l")
plot(fit[[2]]$tau.beta, type="l")
plot(fit[[2]]$rho, type="l")
plot(fit[[2]]$nu, type="l")
plot(fit[[2]]$phi.z, type="l")
plot(fit[[2]]$phi.tau, type="l")
plot(fit[[2]]$phi.w, type="l")
plot(fit[[2]]$beta[, 1], type="l")
plot(fit[[2]]$beta[, 2], type="l")


load("us-all-2.RData")
par(mfrow=c(3, 5))
days1 <- seq(1, 30, by=2)
for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[1]]$tau[, j], type="l", xlab=xlab, main=main)
}
for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[1]]$z[, j], type="l", xlab=xlab, main=main)
}
plot(fit[[1]]$alpha, type="l")
plot(fit[[1]]$tau.alpha, type="l")
plot(fit[[1]]$tau.beta, type="l")
plot(fit[[1]]$rho, type="l")
plot(fit[[1]]$nu, type="l")
plot(fit[[1]]$beta[, 1], type="l")
plot(fit[[1]]$beta[, 2], type="l")

for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[2]]$tau[, j], type="l", xlab=xlab, main=main)
}
for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[2]]$z[, j], type="l", xlab=xlab, main=main)
}

plot(fit[[2]]$alpha, type="l")
plot(fit[[2]]$tau.alpha, type="l")
plot(fit[[2]]$tau.beta, type="l")
plot(fit[[2]]$rho, type="l")
plot(fit[[2]]$nu, type="l")
plot(fit[[2]]$beta[, 1], type="l")
plot(fit[[2]]$beta[, 2], type="l")

load("us-all-27.RData")
par(mfrow=c(3, 5))
days1 <- seq(1, 30, by=2)
for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[1]]$tau[, j], type="l", xlab=xlab, main=main)
}
for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[1]]$z[, j], type="l", xlab=xlab, main=main)
}
plot(fit[[1]]$alpha, type="l")
plot(fit[[1]]$tau.alpha, type="l")
plot(fit[[1]]$tau.beta, type="l")
plot(fit[[1]]$rho, type="l")
plot(fit[[1]]$nu, type="l")
plot(fit[[1]]$phi.z, type="l")
plot(fit[[1]]$phi.tau, type="l")
plot(fit[[1]]$beta[, 1], type="l")
plot(fit[[1]]$beta[, 2], type="l")

for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[2]]$tau[, j], type="l", xlab=xlab, main=main)
}
for (j in days1) {
  xlab <- print(paste("knot 1"))
  main <- print(paste("day", j))
  plot(fit[[2]]$z[, j], type="l", xlab=xlab, main=main)
}

plot(fit[[2]]$alpha, type="l")
plot(fit[[2]]$tau.alpha, type="l")
plot(fit[[2]]$tau.beta, type="l")
plot(fit[[2]]$rho, type="l")
plot(fit[[2]]$nu, type="l")
plot(fit[[2]]$phi.z, type="l")
plot(fit[[2]]$phi.tau, type="l")
plot(fit[[2]]$beta[, 1], type="l")
plot(fit[[2]]$beta[, 2], type="l")

# # X11()
# # boxplot(log(fit.schlather$gpdre))
# # lines(log(r^2))

# # X11()
# # lo<-apply(fit.schlather$yp,1:2,quantile,0.025)
# # me<-apply(fit.schlather$yp,1:2,quantile,0.50)
# # hi<-apply(fit.schlather$yp,1:2,quantile,0.975)

# # plot(y.full.p,me)
# # abline(0,1)
# # print(mean(y.full.p>lo & y.full.p<hi))