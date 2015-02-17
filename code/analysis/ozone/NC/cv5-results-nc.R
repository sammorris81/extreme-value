load("cv-setup.RData")
source("../../R/auxfunctions.R")

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

quant.score.gau <- matrix(NA, length(probs), nsets)
quant.score.t10 <- matrix(NA, length(probs), nsets)
quant.score.t19 <- matrix(NA, length(probs), nsets)
quant.score.t50 <- matrix(NA, length(probs), nsets)
quant.score.t59 <- matrix(NA, length(probs), nsets)

brier.score.gau <- matrix(NA, length(thresholds), nsets)
brier.score.t10 <- matrix(NA, length(thresholds), nsets)
brier.score.t19 <- matrix(NA, length(thresholds), nsets)
brier.score.t50 <- matrix(NA, length(thresholds), nsets)
brier.score.t59 <- matrix(NA, length(thresholds), nsets)

usable <- (25000+1):30000
load("cv5-1NC.RData")
for (d in 1:nsets) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- y[val.idx, ]
  pred.d <- fit.d$yp[usable, , ]
  quant.score.gau[, d] <- QuantScore(pred.d, probs, validate)
  brier.score.gau[, d] <- BrierScore(pred.d, thresholds, validate)
}

load("cv5-2NC.RData")
for (d in 1:nsets) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- y[val.idx, ]
  pred.d <- fit.d$yp[usable, , ]
  quant.score.t10[, d] <- QuantScore(pred.d, probs, validate)
  brier.score.t10[, d] <- BrierScore(pred.d, thresholds, validate)
}

load("cv5-3NC.RData")
for (d in 1:nsets) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- y[val.idx, ]
  pred.d <- fit.d$yp[usable, , ]
  quant.score.t50[, d] <- QuantScore(pred.d, probs, validate)
  brier.score.t50[, d] <- BrierScore(pred.d, thresholds, validate)
}

load("cv5-4NC.RData")
for (d in 1:nsets) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- y[val.idx, ]
  pred.d <- fit.d$yp[usable, , ]
  quant.score.t19[, d] <- QuantScore(pred.d, probs, validate)
  brier.score.t19[, d] <- BrierScore(pred.d, thresholds, validate)
}

load("cv5-5NC.RData")
for (d in 1:nsets) {
  fit.d <- fit[[d]]
  val.idx <- cv.lst[[d]]
  validate <- y[val.idx, ]
  pred.d <- fit.d$yp[usable, , ]
  quant.score.t59[, d] <- QuantScore(pred.d, probs, validate)
  brier.score.t59[, d] <- BrierScore(pred.d, thresholds, validate)
}

savelist <- list(quant.score.gau, quant.score.t10, quant.score.t19, 
             quant.score.t50, quant.score.t59,
		     brier.score.gau, brier.score.t10, brier.score.t19, 
		     brier.score.t50, brier.score.t59,
		     probs, thresholds)

save(savelist, file="cv-scores-nc.RData")

rm(list=ls())
load("cv-setup-nc.RData")
source("../../R/auxfunctions.R")
load("cv-scores-nc.RData")

quant.score.gau <- savelist[[1]]
quant.score.t10 <- savelist[[2]]
quant.score.t19 <- savelist[[3]]
quant.score.t50 <- savelist[[4]]
quant.score.t59 <- savelist[[5]]
brier.score.gau <- savelist[[6]]
brier.score.t10 <- savelist[[7]]
brier.score.t19 <- savelist[[8]]
brier.score.t50 <- savelist[[9]]
brier.score.t59 <- savelist[[10]]
probs <- savelist[[11]]
thresholds <- savelist[[12]]

quant.score.mean <- matrix(NA, 5, length(probs))
brier.score.mean <- matrix(NA, 5, length(thresholds))

quant.score.mean[1, ] <- apply(quant.score.gau, 1, mean)
quant.score.mean[2, ] <- apply(quant.score.t10, 1, mean)
quant.score.mean[3, ] <- apply(quant.score.t19, 1, mean)
quant.score.mean[4, ] <- apply(quant.score.t50, 1, mean)
quant.score.mean[5, ] <- apply(quant.score.t59, 1, mean)

brier.score.mean[1, ] <- apply(brier.score.gau, 1, mean)
brier.score.mean[2, ] <- apply(brier.score.t10, 1, mean)
brier.score.mean[3, ] <- apply(brier.score.t19, 1, mean)
brier.score.mean[4, ] <- apply(brier.score.t50, 1, mean)
brier.score.mean[5, ] <- apply(brier.score.t59, 1, mean)

plot(probs, quant.score.mean[1, ], lty=1, type="b", ylim=c(min(quant.score.mean), max(quant.score.mean)), main="Quantile Scores for ozone analysis (NC only)", xlab="quantile", ylab="score")
for (i in 2:5) {
  lines(probs, quant.score.mean[i, ], lty=i)
  points(probs, quant.score.mean[i, ], pch=i)
}
legend("bottomleft", lty=1:5, pch=1:5, legend=c("Gaussian", "t-1 (T=0.0)", "t-1 (T=0.9)", "t-5 (T=0.0)", "t-5 (T=0.9)"))

# quant.score.mean <- array(NA, dim=c(3, length(probs), 3))
# brier.score.mean <- array(NA, dim=c(3, length(probs), 3))
# for(setting in 1:3){
	# quant.score.mean[1,,setting] <- rowMeans(quant.score.gpd[,,setting])
	# quant.score.mean[2,,setting] <- rowMeans(quant.score.gam[,,setting])
	# quant.score.mean[3,,setting] <- rowMeans(quant.score.mvn[,,setting])
	# brier.score.mean[1,,setting] <- rowMeans(brier.score.gpd[,,setting])
	# brier.score.mean[2,,setting] <- rowMeans(brier.score.gam[,,setting])
	# brier.score.mean[3,,setting] <- rowMeans(brier.score.mvn[,,setting])
# }

# usable <- (burn+1):iters
# lomehi <- c(0.025, 0.5, 0.975)
# #### Setting 1: Threshold: 0.80
# betas.gpd.1 <- apply(betas.gpd[usable,,,1], 2, quantile, probs=lomehi)
# betas.gam.1 <- apply(betas.gam[usable,,,1], 2, quantile, probs=lomehi)
# betas.mvn.1 <- apply(betas.mvn[usable,,,1], 2, quantile, probs=lomehi)

# #### Setting 2: Threshold: 0.95
# betas.gpd.2 <- apply(betas.gpd[usable,,,2], 2, quantile, probs=lomehi)
# betas.gam.2 <- apply(betas.gam[usable,,,2], 2, quantile, probs=lomehi)
# betas.mvn.2 <- apply(betas.mvn[usable,,,2], 2, quantile, probs=lomehi)

# #### Setting 3: Threshold: 0.99
# betas.gpd.3 <- apply(betas.gpd[usable,,,3], 2, quantile, probs=lomehi)
# betas.gam.3 <- apply(betas.gam[usable,,,3], 2, quantile, probs=lomehi)
# betas.mvn.3 <- apply(betas.mvn[usable,,,3], 2, quantile, probs=lomehi)