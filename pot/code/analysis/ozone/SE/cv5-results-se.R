load("cv-setup-se.RData")
source("../../../R/auxfunctions.R")

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

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999)
thresholds <- quantile(y, probs=probs, na.rm=T)
nsets <- 5 # Number of cv sets
nbetas <- 4 # number of betas

quant.score <- array(NA, dim=c(length(probs), nsets, 16))
brier.score <- array(NA, dim=c(length(thresholds), nsets, 16))

beta.0 <- array(NA, dim=c(5000, nsets, 16))
beta.1 <- array(NA, dim=c(5000, nsets, 16))
beta.2 <- array(NA, dim=c(5000, nsets, 16))
beta.3 <- array(NA, dim=c(5000, nsets, 8)) 

# usable <- (25000+1):30000

for (i in 1:16) {
  file <- paste("cv5-", i, "SE.RData", sep="")
  load(file)
  for (d in 1:nsets) {
    fit.d <- fit[[d]]
    val.idx <- cv.lst[[d]]
    validate <- y[val.idx, ]
    pred.d <- fit.d$yp[, , ]
    quant.score[, d, i] <- QuantScore(pred.d, probs, validate)
    brier.score[, d, i] <- BrierScore(pred.d, thresholds, validate)
    beta.0[, d, i] <- fit.d$beta[, 1]
    beta.1[, d, i] <- fit.d$beta[, 2]
    beta.2[, d, i] <- fit.d$beta[, 3]
    if (i <= 8) {
      beta.3[, d, i] <- fit.d$beta[, 4]
    }
  }
}

savelist <- list(quant.score, brier.score,
                 beta.0, beta.1, beta.2, beta.3,
                 probs, thresholds)

save(savelist, file="cv-scores-se.RData")

rm(list=ls())
load("cv-setup-se.RData")
source("../../../R/auxfunctions.R")
load("cv-scores-se.RData")

quant.score <- savelist[[1]]
brier.score <- savelist[[2]]
beta.0 <- savelist[[3]]
beta.1 <- savelist[[4]]
beta.2 <- savelist[[5]]
beta.3 <- savelist[[6]]
probs <- savelist[[7]]
thresholds <- savelist[[8]]

quant.score.mean <- matrix(NA, 16, length(probs))
brier.score.mean <- matrix(NA, 16, length(thresholds))

for (i in 1:16) {
  quant.score.mean[i, ] <- apply(quant.score[, , i], 1, mean)
  brier.score.mean[i, ] <- apply(brier.score[, , i], 1, mean)
}

# plot(probs, quant.score.mean[1, ], lty=1, type="b", ylim=c(min(quant.score.mean), max(quant.score.mean)), main="Quantile Scores for ozone analysis (NC only)", xlab="quantile", ylab="score")
# for (i in 2:5) {
  # lines(probs, quant.score.mean[i, ], lty=i)
  # points(probs, quant.score.mean[i, ], pch=i)
# }
# legend("bottomleft", lty=1:5, pch=1:5, legend=c("Gaussian", "t-1 (T=0.0)", "t-1 (T=0.9)", "t-5 (T=0.0)", "t-5 (T=0.9)"))
