rm(list = ls())
library(compiler)
enableJIT(3)
load("us-all-setup.RData")
source("../../../R/auxfunctions.R")
settings <- read.csv("settings.csv")

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
thresholds <- quantile(Y, probs = probs, na.rm = T)
nsets <- 2 # Number of cv sets
nbetas <- 2 # number of betas

quant.score <- array(NA, dim = c(length(probs), nsets, 74))
brier.score <- array(NA, dim = c(length(thresholds), nsets, 74))

beta.0 <- array(NA, dim = c(5000, nsets, 74))
beta.1 <- array(NA, dim = c(5000, nsets, 74))

phi.z <- array(NA, dim = c(5000, nsets, 74))
phi.w <- array(NA, dim = c(5000, nsets, 74))
phi.tau <- array(NA, dim = c(5000, nsets, 74))

# load("us-all-results.RData")
done <- c(1:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)
# quant.score <- savelist[[1]]
# brier.score <- savelist[[2]]
# beta.0 <- savelist[[3]]
# beta.1 <- savelist[[4]]
# probs <- savelist[[5]]
# thresholds <- savelist[[6]]

for (i in 1:74) {
  file <- paste("results/us-all-", i, ".RData", sep = "")
  cat("start file", file, "\n")
  if (i %in% done) {
    load(file)
    for (d in 1:nsets) {
      if (i == 2) {
        trans <- TRUE
      } else {
        trans <- FALSE
      }
      fit.d <- fit[[d]]
      val.idx <- cv.lst[[d]]
      validate <- Y[val.idx, ]
      pred.d <- fit.d$yp[, , ]
      quant.score[, d, i] <- QuantScore(pred.d, probs, validate, trans = trans)
      brier.score[, d, i] <- BrierScore(pred.d, thresholds, validate, 
                                        trans = trans)
      if (i != 2) {
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

save(savelist, file = "us-all-results.RData")

rm(list = ls())
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

quant.score.mean <- matrix(NA, 74, length(probs))
brier.score.mean <- matrix(NA, 74, length(thresholds))

quant.score.se <- matrix(NA, 74, length(probs))
brier.score.se <- matrix(NA, 74, length(thresholds))

done <- c(1:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)
for (i in 1:74) {
  if (i %in% done) {
    quant.score.mean[i, ] <- apply(quant.score[, , i], 1, mean, na.rm = T)
    quant.score.se[i, ] <- apply(quant.score[, , i], 1, sd) / sqrt(2)
    brier.score.mean[i, ] <- apply(brier.score[, , i], 1, mean, na.rm = T)
    brier.score.se[i, ] <- apply(brier.score[, , i], 1, sd) / sqrt(2)
  }
}

quant.score.mean[c(1:10, 13, 14, 17:19, 21, 23, 25, 26), ]
for (i in 1:length(thresholds)) {
 print(which(quant.score.mean[, i] == min(quant.score.mean[, i], na.rm = T)))
}

for (i in 1:length(thresholds)) {
  print(which(brier.score.mean[, i] == min(brier.score.mean[, i], na.rm = T)))
}

bs.mean.ref.gau <- matrix(NA, nrow = 73, ncol = 11)
qs.mean.ref.gau <- matrix(NA, nrow = 73, ncol = 11)
for (i in 1:73) {
  bs.mean.ref.gau[i, ] <- brier.score.mean[(i + 1), ] / brier.score.mean[1, ]
  qs.mean.ref.gau[i, ] <- quant.score.mean[(i + 1), ] / quant.score.mean[1, ]
}

# find Brier score for methods 36 and 16 (the best performing around 75ppb)
brier.score.site <- matrix(NA, 800, 3)
load('us-all-1.RData')
for (d in 1:2) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- Y[val.idx, ]
  pred.d <- fit.d$yp[, , ]
  brier.score.site[val.idx, 1] <- BrierScoreSite(pred.d, 75, validate)
}

load('us-all-16.RData')
for (d in 1:2) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- Y[val.idx, ]
  pred.d <- fit.d$yp[, , ]
  brier.score.site[val.idx, 2] <- BrierScoreSite(pred.d, 75, validate)
}

load('us-all-36.RData')
for (d in 1:2) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- Y[val.idx, ]
  pred.d <- fit.d$yp[, , ]
  brier.score.site[val.idx, 3] <- BrierScoreSite(pred.d, 75, validate)
}

# get tau such that q(tau) = 75 for each site
library(np)
Y < 75
ozone.quant.site <- apply(Y <= 75, 1, mean, na.rm = T)
sum(ozone.quant.site < 0.5)
sum(ozone.quant.site < 0.7)
sum(ozone.quant.site < 0.76)
sum(ozone.quant.site < 0.78)
sum(ozone.quant.site < 0.82)
sum(ozone.quant.site < 0.83)
quants <- c(0.0, 0.6, 0.8, 0.84, 0.89, 0.91, 0.94, 0.97, 1.01)
in.bin <- rep(0, length(quants) - 1)
sum.bs <- rep(0, length(quants) - 1)
for (q in 1:(length(quants) - 1)) {
  these.q <- which(ozone.quant.site >= quants[q] & 
                   ozone.quant.site < quants[q + 1])
  in.bin[q] <- in.bin[q] + length(these.q)
  sum.bs[q] <- sum.bs[q] + sum(brier.score.site[these.q, 1])
}
mean.bs <- sum.bs / in.bin

these <- which(ozone.quant.site < 1)  # compare when ozone > 75 at least once
par(mfrow=c(1, 2))
yplot <- brier.score.site[, 2] / brier.score.site[, 1]
plot(ozone.quant.site[these], yplot[these], main = "16",
     ylab = "Relative Brier score", 
     xlab = bquote(paste("Marginal ", q^{-1},"(75 ppb)")))
yplot <- brier.score.site[, 3] / brier.score.site[, 1]
plot(ozone.quant.site[these], yplot[these], main = "32",
     ylab = "Relative Brier score", 
     xlab = bquote(paste("Marginal ", q^{-1},"(75 ppb)")))


# get the plotting order for line
plot.ord <- order(ozone.quant.site)

fit.np <- npreg(brier.score.site[, 1] ~ ozone.quant.site, regtype = "ll", 
                bwmethod = "cv.aic", bwtype = "adaptive_nn")
lines(ozone.quant.site[plot.ord], fitted(fit.np)[plot.ord], lty = 2)

# create storage for all 800 sites

# find top two for selected quantiles
score.compare <- bs.mean.ref.gau[, c(1, 6, 9:11)]
idx <- which(score.compare[, 1] == min(score.compare[, 1], na.rm = T))
score.compare[idx, 1]  # 33
idx <- which(score.compare[, 2] == min(score.compare[, 2], na.rm = T))
score.compare[idx, 2]
idx <- which(score.compare[, 3] == min(score.compare[, 3], na.rm = T))
score.compare[idx, 3]
idx <- which(score.compare[, 4] == min(score.compare[, 4], na.rm = T))
score.compare[idx, 4]
idx <- which(score.compare[, 5] == min(score.compare[, 5], na.rm = T))
score.compare[idx, 5]

idx <- order(score.compare[, 1])[2]
score.compare[idx, 1]
idx <- order(score.compare[, 2])[2]
score.compare[idx, 2]
idx <- order(score.compare[, 3])[2]
score.compare[idx, 3]
idx <- order(score.compare[, 4])[2]
score.compare[idx, 4]
idx <- order(score.compare[, 5])[2]
score.compare[idx, 5]
# three main plots (keep max-stable in all for now)
#   time series vs no time-series
#   3 different threshold levels
#   1 knot vs 5 - 10 knots vs 15 knots

these <- 6:11
x.plot <- probs[these]
bg <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1", "gray80")
col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4", "gray16")

# first plot: time series and threshold level (6 lines)
# T = 0
thresh.0.nts  <- c(3, 7, 33, 34, 35, 36, 11, 15) - 1
thresh.0.ts   <- c(51, 54, 57, 60, 63, 66, 69, 72) - 1
bs.mean.ref.gau[thresh.0.nts, ]
bs.mean.ref.gau[thresh.0.ts, ]

# T = 50
thresh.50.nts <- c(4, 8, 38, 39, 40, 41, 12, 16) - 1
thresh.50.ts  <- c(52, 55, 58, 61, 64, 67, 70, 73) - 1
bs.mean.ref.gau[thresh.50.nts, ]
bs.mean.ref.gau[thresh.50.ts, ]

# T = 75
thresh.75.nts <- c(5, 9, 43, 44, 45, 46, 13, 17) - 1
thresh.75.ts  <- c(53, 56, 59, 62, 65, 68, 71, 74) - 1
bs.mean.ref.gau[thresh.75.nts, ]
bs.mean.ref.gau[thresh.75.ts, ]

y.plot <- vector(mode = "list", length=6)
y.plot[[1]] <- apply(bs.mean.ref.gau[thresh.0.nts, these], 2, mean)
y.plot[[2]] <- apply(bs.mean.ref.gau[thresh.50.nts, these], 2, mean)
y.plot[[3]] <- apply(bs.mean.ref.gau[thresh.75.nts, these], 2, mean)
y.plot[[4]] <- apply(bs.mean.ref.gau[thresh.0.ts, these], 2, mean)
y.plot[[5]] <- apply(bs.mean.ref.gau[thresh.50.ts, these], 2, mean)
y.plot[[6]] <- apply(bs.mean.ref.gau[thresh.75.ts, these], 2, mean)

bg <- c("firebrick1", "firebrick1", "firebrick1",
        "dodgerblue1", "dodgerblue1", "dodgerblue1")
col <- c("firebrick4", "firebrick4", "firebrick4",
         "dodgerblue4", "dodgerblue4", "dodgerblue4")
pch <- c(21, 22, 23, 21, 22, 23)
lty <- c(1, 2, 3, 1, 2, 3)
legend <- c("T=0, No Time Series", "T=50, No Time Series",
           "T=75, No Time Series", "T=0, Time Series", "T=50, Time Series",
           "T=75, Time Series")

plot(x.plot, y.plot[[1]], type = "b", lty = 1, ylim = c(0.92, 1),
     bg = bg[1], col = col[1], pch = pch[1], 
     main = "Time series with thresholding",
     ylab = "Relative Brier Score", xlab = "Threshold quantile")
for(i in 2:6) {
  lines(x.plot, y.plot[[i]], type = "b", bg = bg[i], col = col[i], pch = pch[i],
        lty = lty[i])
}
legend("bottomleft", legend = legend, col = col, pch = pch, pt.bg = bg, 
       cex = 1.0, lty = lty, box.lty = 1)

# second plot: time series and number of knots (6 lines)
knots.low.nts  <- c(3, 4, 5) - 1
knots.low.ts   <- c(51, 52, 53) - 1
knots.mid.nts  <- c(7, 8, 9, 33, 38, 43, 34, 39, 44, 35, 40, 45, 36, 41, 46,
                    11, 12, 13) - 1
knots.mid.ts   <- c(54:71) - 1
knots.high.nts <- c(15:17) - 1
knots.high.ts  <- c(72:74) - 1

y.plot <- vector(mode="list", length = 6)
y.plot[[1]] <- apply(bs.mean.ref.gau[knots.low.nts, these], 2, mean)
y.plot[[2]] <- apply(bs.mean.ref.gau[knots.mid.nts, these], 2, mean)
y.plot[[3]] <- apply(bs.mean.ref.gau[knots.high.nts, these], 2, mean)
y.plot[[4]] <- apply(bs.mean.ref.gau[knots.low.ts, these], 2, mean)
y.plot[[5]] <- apply(bs.mean.ref.gau[knots.mid.ts, these], 2, mean)
y.plot[[6]] <- apply(bs.mean.ref.gau[knots.high.ts, these], 2, mean)

bg <- c("firebrick1", "firebrick1", "firebrick1",
        "dodgerblue1", "dodgerblue1", "dodgerblue1")
col <- c("firebrick4", "firebrick4", "firebrick4",
         "dodgerblue4", "dodgerblue4", "dodgerblue4")
pch <- c(21, 22, 23, 21, 22, 23)
lty <- c(1, 2, 3, 1, 2, 3)
legend <- c("K=1, No Time Series", "K=5-10, No Time Series",
           "K=15, No Time Series", "K=1, Time Series", "K=5-10, Time Series",
           "K=15, Time Series")

plot(x.plot, y.plot[[1]], type = "b", lty = 1, ylim = c(0.92, 1),
     bg = bg[1], col = col[1], pch = pch[1], main = "Time series with knots",
     ylab = "Relative Brier Score", xlab = "Threshold quantile")
for(i in 2:6) {
  lines(x.plot, y.plot[[i]], type = "b", bg = bg[i], col = col[i], pch = pch[i],
        lty = lty[i])
}
legend("bottomleft", legend = legend, col = col, pch = pch, pt.bg = bg, 
       cex = 1.0, lty = lty, box.lty = 1)

# third plot: threshold level and number of knots (9 lines)
knots.low.0   <- c(3, 51) - 1
knots.low.50  <- c(4, 52) - 1
knots.low.75  <- c(5, 53) - 1
knots.mid.0   <- c(7, 33:36, 11, 54, 57, 60, 63, 66, 69) - 1
knots.mid.50  <- c(8, 38:41, 12, 55, 58, 61, 64, 67, 70) - 1
knots.mid.75  <- c(9, 43:46, 13, 56, 59, 62, 65, 68, 71) - 1
knots.high.0  <- c(15, 72) - 1
knots.high.50 <- c(16, 73) - 1
knots.high.75 <- c(17, 74) - 1

y.plot <- vector(mode = "list", length=9)
y.plot[[1]] <- apply(bs.mean.ref.gau[knots.low.0, these], 2, mean)
y.plot[[2]] <- apply(bs.mean.ref.gau[knots.low.50, these], 2, mean)
y.plot[[3]] <- apply(bs.mean.ref.gau[knots.low.75, these], 2, mean)
y.plot[[4]] <- apply(bs.mean.ref.gau[knots.mid.0, these], 2, mean)
y.plot[[5]] <- apply(bs.mean.ref.gau[knots.mid.50, these], 2, mean)
y.plot[[6]] <- apply(bs.mean.ref.gau[knots.mid.75, these], 2, mean)
y.plot[[7]] <- apply(bs.mean.ref.gau[knots.high.0, these], 2, mean)
y.plot[[8]] <- apply(bs.mean.ref.gau[knots.high.50, these], 2, mean)
y.plot[[9]] <- apply(bs.mean.ref.gau[knots.high.75, these], 2, mean)

bg <- c("firebrick1", "firebrick1", "firebrick1",
        "dodgerblue1", "dodgerblue1", "dodgerblue1",
        "darkolivegreen1", "darkolivegreen1", "darkolivegreen1")
col <- c("firebrick4", "firebrick4", "firebrick4",
         "dodgerblue4", "dodgerblue4", "dodgerblue4",
         "darkolivegreen4", "darkolivegreen4", "darkolivegreen4")
pch <- c(21, 22, 23, 21, 22, 23, 21, 22, 23)
lty <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
legend <- c("K=1, T=0", "K=1, T=50", "K=1, T=75",
            "K=5-10, T=0", "K=5-10, T=50", "K=5-10, T=75",
            "K=15, T=0", "K=15, T=50", "K=15, T=75")

plot(x.plot, y.plot[[1]], type = "b", lty = 1, ylim = c(0.92, 1),
     bg = bg[1], col = col[1], pch = pch[1], main = "Knots and thresholding",
     ylab = "Relative Brier Score", xlab = "Threshold quantile")
for(i in 2:9) {
  lines(x.plot, y.plot[[i]], type = "b", bg = bg[i], col = col[i], pch = pch[i],
        lty = lty[i])
}
legend("bottomleft", legend = legend, col = col, pch = pch, pt.bg = bg, 
       cex = 1.0, lty = lty, box.lty = 1)

bg <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1", "gray80")
col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4", "gray16")

# another set of plots, 1 time series, 1 no time series
# bs-ozone.pdf
these <- 6:11
x.plot <- probs[these]
bg <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1", "gray80")
col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4", "gray16")

# time series
y.plot <- vector(mode = "list", length = 9)
y.plot[[1]] <- bs.mean.ref.gau[50, these]
y.plot[[2]] <- bs.mean.ref.gau[51, these]
y.plot[[3]] <- bs.mean.ref.gau[52, these]
y.plot[[4]] <- bs.mean.ref.gau[60, these]
y.plot[[5]] <- bs.mean.ref.gau[61, these]
y.plot[[6]] <- bs.mean.ref.gau[62, these]
y.plot[[7]] <- bs.mean.ref.gau[71, these]
y.plot[[8]] <- bs.mean.ref.gau[72, these]
y.plot[[9]] <- bs.mean.ref.gau[73, these]

bg <- c("firebrick1", "firebrick1", "firebrick1",
        "dodgerblue1", "dodgerblue1", "dodgerblue1",
        "darkolivegreen1", "darkolivegreen1", "darkolivegreen1")
col <- c("firebrick4", "firebrick4", "firebrick4",
         "dodgerblue4", "dodgerblue4", "dodgerblue4",
         "darkolivegreen4", "darkolivegreen4", "darkolivegreen4")
pch <- c(21, 22, 23, 21, 22, 23, 21, 22, 23)
lty <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)

# Panel for paper
quartz(width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5.1, 5.1, 4.1, 2.1))

plot(x.plot, y.plot[[1]], type="b", lty = 1, ylim = c(0.9, 1),
     bg = bg[1], col = col[1], pch = pch[1], # main = "Time Series Models",
     ylab = "Relative Brier Score", xlab = "Threshold quantile",
     cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main=2)
for(i in 2:9) {
  lines(x.plot, y.plot[[i]], type = "b", bg = bg[i], col = col[i], pch = pch[i],
        lty = lty[i], cex = 1.5)
}
# abline(h=1, lty=2)
legend <- c("K=1, T=0", "K=1, T=50", "K=1, T=75",
            "K=7, T=0", "K=7, T=50", "K=7, T=75",
            "K=15, T=0", "K=15, T=50", "K=15, T=75")
legend("bottomleft", legend = legend, col = col, pch = pch, pt.bg = bg, 
       cex = 1.0, lty = lty, box.lty = 1)

# non time-series
y.plot <- vector(mode = "list", length = 9)
y.plot[[1]] <- bs.mean.ref.gau[2, these]
y.plot[[2]] <- bs.mean.ref.gau[3, these]
y.plot[[3]] <- bs.mean.ref.gau[4, these]
y.plot[[4]] <- bs.mean.ref.gau[33, these]
y.plot[[5]] <- bs.mean.ref.gau[38, these]
y.plot[[6]] <- bs.mean.ref.gau[43, these]
y.plot[[7]] <- bs.mean.ref.gau[14, these]
y.plot[[8]] <- bs.mean.ref.gau[15, these]
y.plot[[9]] <- bs.mean.ref.gau[16, these]

bg <- c("firebrick1", "firebrick1", "firebrick1",
        "dodgerblue1", "dodgerblue1", "dodgerblue1",
        "darkolivegreen1", "darkolivegreen1", "darkolivegreen1")
col <- c("firebrick4", "firebrick4", "firebrick4",
         "dodgerblue4", "dodgerblue4", "dodgerblue4",
         "darkolivegreen4", "darkolivegreen4", "darkolivegreen4")
pch <- c(21, 22, 23, 21, 22, 23, 21, 22, 23)
lty <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)

plot(x.plot, y.plot[[1]], type="b", lty=1, ylim=c(0.9, 1),
     bg=bg[1], col=col[1], pch=pch[1], # main="Non Time Series Models",
     ylab="Relative Brier Score", xlab="Threshold quantile",
     cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=2)
for(i in 2:9) {
  lines(x.plot, y.plot[[i]], type="b", bg=bg[i], col=col[i], pch=pch[i],
        lty=lty[i], cex=1.5)
}
# abline(h=1, lty=2)
# legend <- c("K=1, T=0", "K=1, T=50", "K=1, T=75",
#             "K=7, T=0", "K=7, T=50", "K=7, T=75",
#             "K=15, T=0", "K=15, T=50", "K=15, T=75")
# legend("bottomleft", legend=legend, col=col, pch=pch, pt.bg=bg, cex=1.0,
#        lty=lty, box.lty=1)


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
bg <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1", "gray80")
col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4", "gray16")
these.probs <- 1:11
xplot <- probs[these.probs]
plot(xplot, bs.mean.ref.gau[1, these.probs], type="b", ylim=c(0.755, 0.95), pch=21,
     col=col[1], bg=bg[1], ylab="Relative Brier score", xlab="Threshold quantile", lty=1,
     # main="Select ozone Brier scores",
     cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex=1.3)
lines(xplot, bs.mean.ref.gau[5, these.probs], type="b", pch=21, col=col[2], bg=bg[2],
      cex=1.3, lty=1)
lines(xplot, bs.mean.ref.gau[9, these.probs], type="b", pch=21, col=col[3], bg=bg[3],
      cex=1.3, lty=1)
lines(xplot, bs.mean.ref.gau[13, these.probs], type="b", pch=21, col=col[4], bg=bg[4],
      cex=1.3, lty=1)
lines(xplot[3:11], bs.mean.ref.gau[25, 3:11], type="b", pch=23, col=col[5], bg=bg[5],
      cex=1.3, lty=1)
abline(h=1, lty=2)
legend("topleft", legend=c("Skew-t, K=1, T=0", "Skew-t, K=5, T=0", "Skew-t, K=10, T=0",
       "Skew-t, K=15, T=0", "Max-stable"),
       pch=c(21, 21, 21, 21, 23), lty=c(1, 1, 1, 1, 1), cex=1.3,
       pt.bg=c(bg[1], bg[2], bg[3], bg[4], bg[5]),
       col=c(col[1], col[2], col[3], col[4], col[5]))

bg <- c("firebrick1", "dodgerblue1", "darkolivegreen1", "orange1", "gray80")
col <- c("firebrick4", "dodgerblue4", "darkolivegreen4", "orange4", "gray16")
these.probs <- 1:11
xplot <- probs[these.probs]
plot(xplot, qs.mean.ref.gau[1, these.probs], type="b", ylim=c(0.5, 1.1), pch=21,
     col=col[1], bg=bg[1], ylab="Relative quantile score", xlab="Threshold quantile", lty=1,
     # main="Select ozone Brier scores",
     cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex=1.3)
lines(xplot, qs.mean.ref.gau[5, these.probs], type="b", pch=21, col=col[2], bg=bg[2],
      cex=1.3, lty=1)
lines(xplot, qs.mean.ref.gau[9, these.probs], type="b", pch=21, col=col[3], bg=bg[3],
      cex=1.3, lty=1)
lines(xplot, qs.mean.ref.gau[13, these.probs], type="b", pch=21, col=col[4], bg=bg[4],
      cex=1.3, lty=1)
lines(xplot[3:11], qs.mean.ref.gau[25, 3:11], type="b", pch=23, col=col[5], bg=bg[5],
      cex=1.3, lty=1)
abline(h=1, lty=2)
legend("bottomleft", legend=c("Skew-t, K=1, T=0", "Skew-t, K=5, T=0", "Skew-t, K=10, T=0",
       "Skew-t, K=15, T=0", "Max-stable"),
       pch=c(21, 21, 21, 21, 23), lty=c(1, 1, 1, 1, 1), cex=1.3,
       pt.bg=c(bg[1], bg[2], bg[3], bg[4], bg[5]),
       col=c(col[1], col[2], col[3], col[4], col[5]))



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



# posterior predictions
rm(list=ls())
load("../ozone_data.RData")

# preprocessing
x <- x / 1000
y <- y / 1000

# get the locations of
S.p <- expand.grid(x, y)
keep.these <- (S.p[, 1] > 1.03 & S.p[, 1] < 1.7) &
  (S.p[, 2] > -0.96 & S.p[, 2] < -0.40)
S.p     <- S.p[keep.these, ]
nx <- length(unique(S.p[, 1]))
ny <- length(unique(S.p[, 2]))
nt <- dim(Y)[2]

# load results
threshold <- 75
load('us-all-pred-1.RData')
yp <- y.pred
np <- dim(yp)[2]
set.1.95 <- apply(yp, c(2), quantile, probs=0.95)
set.1.99 <- apply(yp, c(2), quantile, probs=0.99)
set.1.p.below <- matrix(0, np, nt)
for (i in 1:np) { for (t in 1:nt) {
  set.1.p.below[i, t] <- mean(yp[, i, t] <= threshold)
} }
set.1.p.0 <- rep(0, np)
for(i in 1:np) {
  set.1.p.0[i] <- prod(set.1.p.below[i, ])
}
set.1.p.1 <- rep(0, np)
for (i in 1:np) { for (t in 1:nt) {
  set.1.p.1[i] <- set.1.p.1[i] + prod(set.1.p.below[i, -t]) *
                  (1 - set.1.p.below[i, t])
} }
set.1.p.2 <- rep(0, np)
for(i in 1:np) { for (t in 1:(nt - 1)) {
  for (s in (t+1):nt) {
    set.1.p.2[i] <- set.1.p.2[i] + prod(set.1.p.below[i, -c(s,t)]) *
                    prod(1 - set.1.p.below[i, c(s, t)])
  }
}}
set.1.p.atleast1 <- 1 - set.1.p.0
set.1.p.atleast2 <- 1 - (set.1.p.0 + set.1.p.1)
set.1.p.atleast3 <- 1 - (set.1.p.0 + set.1.p.1 + set.1.p.2)


# 1 knot - No Time Series - T = 0
load('us-all-pred-3.RData')
yp <- y.pred
np <- dim(yp)[2]
set.3.95 <- apply(yp, c(2), quantile, probs=0.95)
set.3.99 <- apply(yp, c(2), quantile, probs=0.99)
set.3.p.below <- matrix(0, np, nt)
for (i in 1:np) { for (t in 1:nt) {
  set.3.p.below[i, t] <- mean(yp[, i, t] <= threshold)
} }
set.3.p.0 <- rep(0, np)
for(i in 1:np) {
  set.3.p.0[i] <- prod(set.3.p.below[i, ])
}
set.3.p.1 <- rep(0, np)
for (i in 1:np) { for (t in 1:nt) {
  set.3.p.1[i] <- set.3.p.1[i] + prod(set.3.p.below[i, -t]) *
                  (1 - set.3.p.below[i, t])
} }
set.3.p.2 <- rep(0, np)
for(i in 1:np) { for (t in 1:(nt - 1)) {
  for (s in (t+1):nt) {
    set.3.p.2[i] <- set.3.p.2[i] + prod(set.3.p.below[i, -c(s,t)]) *
                    prod(1 - set.3.p.below[i, c(s, t)])
  }
}}
set.3.p.atleast1 <- 1 - set.3.p.0
set.3.p.atleast2 <- 1 - (set.3.p.0 + set.3.p.1)
set.3.p.atleast3 <- 1 - (set.3.p.0 + set.3.p.1 + set.3.p.2)



# Skew-t - No Time series - T = 50
load('us-all-pred-8.RData')
yp <- y.pred
np <- dim(yp)[2]
set.8.95 <- apply(yp, c(2), quantile, probs=0.95)
set.8.99 <- apply(yp, c(2), quantile, probs=0.99)
set.8.p.below <- matrix(0, np, nt)
for (i in 1:np) { for (t in 1:nt) {
  set.8.p.below[i, t] <- mean(yp[, i, t] <= threshold)
} }
set.8.p.0 <- rep(0, np)
for(i in 1:np) {
  set.8.p.0[i] <- prod(set.8.p.below[i, ])
}
set.8.p.1 <- rep(0, np)
for (i in 1:np) { for (t in 1:nt) {
  set.8.p.1[i] <- set.8.p.1[i] + prod(set.8.p.below[i, -t]) *
    (1 - set.8.p.below[i, t])
} }
set.8.p.2 <- rep(0, np)
for(i in 1:np) { for (t in 1:(nt - 1)) {
  for (s in (t+1):nt) {
    set.8.p.2[i] <- set.8.p.2[i] + prod(set.8.p.below[i, -c(s,t)]) *
      prod(1 - set.8.p.below[i, c(s, t)])
  }
}}
set.8.p.atleast1 <- 1 - set.8.p.0
set.8.p.atleast2 <- 1 - (set.8.p.0 + set.8.p.1)
set.8.p.atleast3 <- 1 - (set.8.p.0 + set.8.p.1 + set.8.p.2)

# 6 knots - Time series - T = 75
load('us-all-pred-59.RData')
yp <- y.pred
np <- dim(yp)[2]
set.59.95 <- apply(yp, c(2), quantile, probs=0.95)
set.59.99 <- apply(yp, c(2), quantile, probs=0.99)
set.59.p.below <- matrix(0, np, nt)
for (i in 1:np) { for (t in 1:nt) {
  set.59.p.below[i, t] <- mean(yp[, i, t] <= threshold)
} }
set.59.p.0 <- rep(0, np)
for(i in 1:np) {
  set.59.p.0[i] <- prod(set.59.p.below[i, ])
}
set.59.p.1 <- rep(0, np)
for (i in 1:np) { for (t in 1:nt) {
  set.59.p.1[i] <- set.59.p.1[i] + prod(set.59.p.below[i, -t]) *
    (1 - set.59.p.below[i, t])
} }
set.59.p.2 <- rep(0, np)
for(i in 1:np) { for (t in 1:(nt - 1)) {
  for (s in (t+1):nt) {
    set.59.p.2[i] <- set.59.p.2[i] + prod(set.59.p.below[i, -c(s,t)]) *
      prod(1 - set.59.p.below[i, c(s, t)])
  }
}}
set.59.p.atleast1 <- 1 - set.59.p.0
set.59.p.atleast2 <- 1 - (set.59.p.0 + set.59.p.1)
set.59.p.atleast3 <- 1 - (set.59.p.0 + set.59.p.1 + set.59.p.2)

# 10 knots - Time series - T = 75
load('us-all-pred-71.RData')
yp <- y.pred
np <- dim(yp)[2]
set.71.95 <- apply(yp, c(2), quantile, probs=0.95)
set.71.99 <- apply(yp, c(2), quantile, probs=0.99)
set.71.p.below <- matrix(0, np, nt)
for (i in 1:np) { for (t in 1:nt) {
  set.71.p.below[i, t] <- mean(yp[, i, t] <= threshold)
} }
set.71.p.0 <- rep(0, np)
for(i in 1:np) {
  set.71.p.0[i] <- prod(set.71.p.below[i, ])
}
set.71.p.1 <- rep(0, np)
for (i in 1:np) { for (t in 1:nt) {
  set.71.p.1[i] <- set.71.p.1[i] + prod(set.71.p.below[i, -t]) *
                   (1 - set.71.p.below[i, t])
} }
set.71.p.2 <- rep(0, np)
for(i in 1:np) { for (t in 1:(nt - 1)) {
  for (s in (t+1):nt) {
    set.71.p.2[i] <- set.71.p.2[i] + prod(set.71.p.below[i, -c(s,t)]) *
                     prod(1 - set.71.p.below[i, c(s, t)])
  }
}}
set.71.p.atleast1 <- 1 - set.71.p.0
set.71.p.atleast2 <- 1 - (set.71.p.0 + set.71.p.1)
set.71.p.atleast3 <- 1 - (set.71.p.0 + set.71.p.1 + set.71.p.2)
rm(y.pred)
rm(yp)
save.image(file="predict-maps.RData")

# make the prediction maps
rm(list=ls())
load("predict-maps.RData")
library(fields)
library(SpatialTools)

# probability at least one day exceeds
quartz(width=11, height=8)
par(mfrow=c(2, 3))
quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.1.p.atleast1), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.3.p.atleast1), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.8.p.atleast1), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.59.p.atleast1), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.71.p.atleast1), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

# probability at least two days exceed
quartz(width=12, height=8)
par(mfrow=c(2, 3))
quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.1.p.atleast2), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(0, 1),
           main="Gaus - No Time Series, T=0")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.3.p.atleast2), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(0, 1),
           main="Skew-t, K=1, No Time Series, T=0")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.8.p.atleast2), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.59.p.atleast2), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.71.p.atleast2), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(0, 1),
           main="Sym-t, K=10, Time Series, T=75")
lines(borders/1000)

# probability at least three days exceed
quartz(width=12, height=8)
par(mfrow=c(2, 3))
quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.1.p.atleast3), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.3.p.atleast3), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.8.p.atleast3), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.59.p.atleast3), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.71.p.atleast3), nx=nx, ny=ny,
           yaxt="n", xaxt="n")
lines(borders/1000)

# 95th quantiles
quartz(width=12, height=8)
par(mfrow=c(2, 3))
quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.1.95), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(55, 100),
           main="Gaus - No Time Series, T=0")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.3.95), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(55, 100),
           main="Skew-t, K=1, No Time Series, T=0")
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.8.95), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(55, 100))
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.59.95), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(55, 100))
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.71.95), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(55, 100),
           main="Sym-t, K=10, Time Series, T=75")
lines(borders/1000)

# 99th quantiles
quartz(width=6, height=9)
par(mfrow=c(3, 2), mar = c(3, 2, 3, 2))
quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.1.99), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(60, 120),
           main="(a) Gaussian", cex=1.5)
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.3.99), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(60, 120),
           main="(b) Skew-t, K=1, T=0, No Time Series, ", cex=1.5)
lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.8.99), nx=nx, ny=ny,
           yaxt="n", xaxt="n",zlim=c(60, 120),
           main="(c) Skew-t, K=5, T=50, No Time Series")
lines(borders/1000)
#
# quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.59.99), nx=nx, ny=ny,
#            yaxt="n", xaxt="n",zlim=c(60, 120),
#            main="Skew-t, K=6, Time Series, T=75")
# lines(borders/1000)

quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(set.71.99), nx=nx, ny=ny,
           yaxt="n", xaxt="n", zlim=c(60, 120),
           main="(d) Sym-t, K=10, T=75, Time Series", cex=1.5)
lines(borders/1000)

diff.71.1 = set.71.99 - set.1.99
quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(diff.71.1), nx=nx, ny=ny,
           yaxt="n", xaxt="n",
           main="(e) Difference of (d) - (a)", cex=1.5)
lines(borders/1000)

diff.71.3 = set.71.99 - set.3.99
quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(diff.71.3), nx=nx, ny=ny,
           yaxt="n", xaxt="n",
           main="(f) Difference of (d) - (b)", cex=1.5)
lines(borders/1000)







# zlim=c(0, 122)
# quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(CMAQ.p[, 5]), zlim=zlim, nx=nx, ny=ny)
# lines(borders/1000)

# posterior maps: 95th quantile

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

# get an idea of two sites that are close to one another vs far apart
rm(list=ls())
library(fields)
load("us-all-setup.RData")
source("../../../R/auxfunctions.R")
settings <- read.csv("settings.csv")
no.na <- which(rowSums(is.na(Y)) == 0)
indices <- 1:nrow(S)

# remove the NAs
S <- S[no.na, ]
indices <- indices[no.na]

# find the sites that are close to one another
d <- rdist(S) 
diag(d) <- 0

close <- which(d < 0.012 & d > 0, arr.ind = TRUE)
close.these <- duplicated(rowSums(close))
close <- close[!close.these, ]
far <- which(d > 4.46, arr.ind = TRUE)
far.these <- duplicated(rowSums(far))
far <- far[!far.these, ]

quant.close.1 <- quant.close.2 <- quant.far.1 <- quant.far.2 <- rep(NA, ncol(Y))
for (i in 1:nrow(close)) {
	idx.1 <- indices[close[i, 1]]
	idx.2 <- indices[close[i, 2]]
	for (j in 1:ncol(Y)) {
		quant.close.1[j] <- mean(Y[idx.1, j] > Y[idx.1, -j])
		quant.close.2[j] <- mean(Y[idx.2, j] > Y[idx.2, -j])
	}
	if (i == 1) {
		plot(quant.close.1, quant.close.2)
	} else {
		points(quant.close.1, quant.close.2)
	}
}

# plot(S)
close <- close[1, ]
indices <- indices[close]


load("us-all-setup.RData")

close <- c(47, 215)
far <- c(698, 215)

plot(S, type = 'n')
lines(borders / 1000)
points(S[close, ], ylim=c(-1.6, 1), xlim=c(-2.2, 2.2))
points(S[far, ], ylim=c(-1.6, 1), xlim=c(-2.2, 2.2))
rdist(S[close, ] * 1000)[1, 2]
rdist(S[far, ] * 1000)[1, 2]

quant.close.1 <- quant.close.2 <- quant.far.1 <- quant.far.2 <- rep(NA, ncol(Y))
for (i in 1:ncol(Y)) {
	quant.close.1[i] <- mean(Y[close[1], i] > Y[close[1], -i])
	quant.close.2[i] <- mean(Y[close[2], i] > Y[close[2], -i])
	quant.far.1[i] <- mean(Y[far[1], i] > Y[far[1], -i])
	quant.far.2[i] <- mean(Y[far[2], i] > Y[far[2], -i])
}

plot(t(Y[close, ]), main = "",
     xlab = paste("Site", close[1]), ylab = paste("Site", close[2]))
plot(t(Y[far, ]), main = "",
     xlab = paste("Site", far[1]), ylab = paste("Site", far[2]))   

quartz(width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5.1, 5.1, 4.1, 2.1))
plot(quant.close.1, quant.close.2, main = "12 km apart", cex.lab = 1.5,
     xlab = "Columbus, OH", ylab = "Columbus, OH")
plot(quant.far.1, quant.far.2, main = "3,143 km apart", cex.lab = 1.5, 
     xlab = "Los Angeles, CA", ylab = "Columbus, OH")

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