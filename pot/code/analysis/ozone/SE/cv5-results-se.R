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

usable <- (25000+1):30000

for (i in 1:16) {
  file <- paste("cv5-", i, "SE.RData", sep="")
  load(file)
  for (d in 1:nsets) {
    fit.d <- fit[[d]]
    val.idx <- cv.lst[[d]]
    validate <- y[val.idx, ]
    pred.d <- fit.d$yp[usable, , ]
    quant.score[, d, i] <- QuantScore(pred.d, probs, validate)
    brier.score[, d, i] <- BrierScore(pred.d, thresholds, validate)
    beta.0[, d, i] <- fit.d$beta[usable, 1]
    beta.1[, d, i] <- fit.d$beta[usable, 2]
    beta.2[, d, i] <- fit.d$beta[usable, 3]
    if (i <= 8) {
      beta.3[, d, i] <- fit.d$beta[usable, 4]
    }
  }
}

savelist <- list(quant.score, brier.score,
                 beta.0, beta.1, beta.2, beta.3,
                 probs, thresholds)

save(savelist, file="cv-scores-se.RData")

# rm(list=ls())
# load("cv-setup-se.RData")
# source("../../R/auxfunctions.R")
# load("cv-scores-se.RData")

# quant.score.gau <- savelist[[1]]
# quant.score.t10 <- savelist[[2]]
# quant.score.t19 <- savelist[[3]]
# quant.score.t50 <- savelist[[4]]
# quant.score.t59 <- savelist[[5]]
# brier.score.gau <- savelist[[6]]
# brier.score.t10 <- savelist[[7]]
# brier.score.t19 <- savelist[[8]]
# brier.score.t50 <- savelist[[9]]
# brier.score.t59 <- savelist[[10]]
# probs <- savelist[[11]]
# thresholds <- savelist[[12]]

# quant.score.mean <- matrix(NA, 5, length(probs))
# brier.score.mean <- matrix(NA, 5, length(thresholds))

# quant.score.mean[1, ] <- apply(quant.score.gau, 1, mean)
# quant.score.mean[2, ] <- apply(quant.score.t10, 1, mean)
# quant.score.mean[3, ] <- apply(quant.score.t19, 1, mean)
# quant.score.mean[4, ] <- apply(quant.score.t50, 1, mean)
# quant.score.mean[5, ] <- apply(quant.score.t59, 1, mean)

# brier.score.mean[1, ] <- apply(brier.score.gau, 1, mean)
# brier.score.mean[2, ] <- apply(brier.score.t10, 1, mean)
# brier.score.mean[3, ] <- apply(brier.score.t19, 1, mean)
# brier.score.mean[4, ] <- apply(brier.score.t50, 1, mean)
# brier.score.mean[5, ] <- apply(brier.score.t59, 1, mean)

# plot(probs, quant.score.mean[1, ], lty=1, type="b", ylim=c(min(quant.score.mean), max(quant.score.mean)), main="Quantile Scores for ozone analysis (NC only)", xlab="quantile", ylab="score")
# for (i in 2:5) {
  # lines(probs, quant.score.mean[i, ], lty=i)
  # points(probs, quant.score.mean[i, ], pch=i)
# }
# legend("bottomleft", lty=1:5, pch=1:5, legend=c("Gaussian", "t-1 (T=0.0)", "t-1 (T=0.9)", "t-5 (T=0.0)", "t-5 (T=0.9)"))
