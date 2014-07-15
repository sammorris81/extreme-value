#########################################################################
# A simulation study to determine if a thresholded method improves
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
# 1 - Gaussian
# 2 - skew-Gaussian (alpha = 5)
# 3 - skew-t (alpha = 5)
# 4 - skew-t w/partition (5 knots)
# 5 - 1/2 Gaussian, 1/2 t
# 6 - 1/2 Gaussian, 1/2 t w/partition (5 knots)
# 7 - 1/2 Gaussian (range = 0.40), 1/2 t (range = 0.10)
# 8 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
#
# analysis methods:
#  1 - Gaussian
#  2 - t
#  3 - t w/partition (5 knots)
#  4 - skew-Gaussian
#  5 - skew-t
#  6 - skew-t w/partition (5 knots)
#  7 - t (thresh = 0.90)
#  8 - t w/partition (5 knots, thresh = 0.90)
#  9 - skew-t (thresh=0.90)
#
#########################################################################

rm(list=ls())
load("simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- 5
nsettings <- dim(y)[4]
nmethods <- 9

source("../../R/auxfunctions.R")	# Included for easy access if we need to change score functions

# results should include
#   - coverage for all parameters
#   - quantile score plots for each data setting
#
# result RData files are numbered by datasetting
# setting a: methods 1 - 5
# setting b: methods 6 - 9
#
# results do not have burnin
# fit.1[[2]] are the results for method: Gaussian on the second dataset

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999)

quant.score <- array(NA, dim=c(length(probs), nsets, nmethods, nsettings))
brier.score <- array(NA, dim=c(length(probs), nsets, nmethods, nsettings))

# storage for the interval endpoints
intervals <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
beta.0 <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
beta.1 <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
beta.2 <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
tau.alpha <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
tau.beta <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
rho <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
nu <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
alpha <- array(NA, dim=c(length(intervals), nsets, nmethods, nsettings))
# not all methods use skew or multiple partitions
z.alpha <- array(NA, dim=c(length(intervals), nsets, 4, nsettings))
avgparts <- array(NA, dim=c(length(intervals), nsets, 3, nsettings))

iters <- 20000; burn <- 10000
for(setting in 1:nsettings){
  dataset <- paste(setting,"-a.RData", sep="")
  load(dataset)
  
  for(d in 1:nsets){  # fit.1 is gaussian, fit.2 is t, etc.
  	thresholds <- quantile(y[, , d, setting], probs=probs, na.rm=T)
    validate <- y.validate[, , d]
    
    pred <- fit.1$yp  # gaussian
    quant.score[, d, 1, setting] <- QuantScore(pred, probs, validate) 
    brier.score[, d, 1, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 1, setting] <- quantile(fit.1$beta[, 1], probs=intervals)
    beta.1[, d, 1, setting] <- quantile(fit.1$beta[, 2], probs=intervals)
    beta.2[, d, 1, setting] <- quantile(fit.1$beta[, 3], probs=intervals)
    tau.alpha[, d, 1, setting] <- quantile(fit.1$tau.alpha, probs=intervals)
    tau.beta[, d, 1, setting] <- quantile(fit.1$tau.beta, probs=intervals)
    rho[, d, 1, setting] <- quantile(fit.1$rho, probs=intervals)
    nu[, d, 1, setting] <- quantile(fit.1$nu, probs=intervals)
    alpha[, d, 1, setting] <- quantile(fit.1$alpha, probs=intervals)
    
    pred <- fit.2$yp  # t1
    quant.score[, d, 2, setting] <- QuantScore(pred, probs, validate)
    brier.score[, d, 2, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 2, setting] <- quantile(fit.2$beta[, 1], probs=intervals)
    beta.1[, d, 2, setting] <- quantile(fit.2$beta[, 2], probs=intervals)
    beta.2[, d, 2, setting] <- quantile(fit.2$beta[, 3], probs=intervals)
    tau.alpha[, d, 2, setting] <- quantile(fit.2$tau.alpha, probs=intervals)
    tau.beta[, d, 2, setting] <- quantile(fit.2$tau.beta, probs=intervals)
    rho[, d, 2, setting] <- quantile(fit.2$rho, probs=intervals)
    nu[, d, 2, setting] <- quantile(fit.2$nu, probs=intervals)
    alpha[, d, 2, setting] <- quantile(fit.2$alpha, probs=intervals)
    
    pred <- fit.3$yp  # t5
    quant.score[, d, 3, setting] <- QuantScore(pred, probs, validate)
    brier.score[, d, 3, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 3, setting] <- quantile(fit.3$beta[, 1], probs=intervals)
    beta.1[, d, 3, setting] <- quantile(fit.3$beta[, 2], probs=intervals)
    beta.2[, d, 3, setting] <- quantile(fit.3$beta[, 3], probs=intervals)
    tau.alpha[, d, 3, setting] <- quantile(fit.3$tau.alpha, probs=intervals)
    tau.beta[, d, 3, setting] <- quantile(fit.3$tau.beta, probs=intervals)
    rho[, d, 3, setting] <- quantile(fit.3$rho, probs=intervals)
    nu[, d, 3, setting] <- quantile(fit.3$nu, probs=intervals)
    alpha[, d, 3, setting] <- quantile(fit.3$alpha, probs=intervals)
    avgparts[, d, 1, setting] <- mean(fit.3$avgparts)
    
    pred <- fit.4$yp  # skew-gaussian
    quant.score[, d, 4, setting] <- QuantScore(pred, probs, validate)
    brier.score[, d, 4, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 4, setting] <- quantile(fit.4$beta[, 1], probs=intervals)
    beta.1[, d, 4, setting] <- quantile(fit.4$beta[, 2], probs=intervals)
    beta.2[, d, 4, setting] <- quantile(fit.4$beta[, 3], probs=intervals)
    tau.alpha[, d, 4, setting] <- quantile(fit.4$tau.alpha, probs=intervals)
    tau.beta[, d, 4, setting] <- quantile(fit.4$tau.beta, probs=intervals)
    rho[, d, 4, setting] <- quantile(fit.4$rho, probs=intervals)
    nu[, d, 4, setting] <- quantile(fit.4$nu, probs=intervals)
    alpha[, d, 4, setting] <- quantile(fit.4$alpha, probs=intervals)
    z.alpha[, d, 1, setting] <- quantile(fit.4$z.alpha, probs=intervals)
    
    pred <- fit.5$yp  # skew-t1
    quant.score[, d, 5, setting] <- QuantScore(pred, probs, validate) 
    brier.score[, d, 5, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 5, setting] <- quantile(fit.5$beta[, 1], probs=intervals)
    beta.1[, d, 5, setting] <- quantile(fit.5$beta[, 2], probs=intervals)
    beta.2[, d, 5, setting] <- quantile(fit.5$beta[, 3], probs=intervals)
    tau.alpha[, d, 5, setting] <- quantile(fit.5$tau.alpha, probs=intervals)
    tau.beta[, d, 5, setting] <- quantile(fit.5$tau.beta, probs=intervals)
    rho[, d, 5, setting] <- quantile(fit.5$rho, probs=intervals)
    nu[, d, 5, setting] <- quantile(fit.5$nu, probs=intervals)
    alpha[, d, 5, setting] <- quantile(fit.5$alpha, probs=intervals)
    z.alpha[, d, 2, setting] <- quantile(fit.5$z.alpha, probs=intervals)
  }
  cat("dataset", dataset, "-a \n")
  dataset <- paste(setting,"-b.RData", sep="")
  for (d in 1:nsets){
  	thresholds <- quantile(y[, , d, setting], probs=probs, na.rm=T)
    validate <- y.validate[, , d]
    
    pred <- fit.1$yp  # skew-t5
    quant.score[, d, 6, setting] <- QuantScore(pred, probs, validate) 
    brier.score[, d, 6, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 6, setting] <- quantile(fit.1$beta[, 1], probs=intervals)
    beta.1[, d, 6, setting] <- quantile(fit.1$beta[, 2], probs=intervals)
    beta.2[, d, 6, setting] <- quantile(fit.1$beta[, 3], probs=intervals)
    tau.alpha[, d, 6, setting] <- quantile(fit.1$tau.alpha, probs=intervals)
    tau.beta[, d, 6, setting] <- quantile(fit.1$tau.beta, probs=intervals)
    rho[, d, 6, setting] <- quantile(fit.1$rho, probs=intervals)
    nu[, d, 6, setting] <- quantile(fit.1$nu, probs=intervals)
    alpha[, d, 6, setting] <- quantile(fit.1$alpha, probs=intervals)
    z.alpha[, d, 3, setting] <- quantile(fit.1$z.alpha, probs=intervals)
    avgparts[, d, 2, setting] <- mean(fit.1$avgparts)
    
    pred <- fit.2$yp  # t1 (T=0.90)
    quant.score[, d, 7, setting] <- QuantScore(pred, probs, validate) 
    brier.score[, d, 7, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 7, setting] <- quantile(fit.2$beta[, 1], probs=intervals)
    beta.1[, d, 7, setting] <- quantile(fit.2$beta[, 2], probs=intervals)
    beta.2[, d, 7, setting] <- quantile(fit.2$beta[, 3], probs=intervals)
    tau.alpha[, d, 7, setting] <- quantile(fit.2$tau.alpha, probs=intervals)
    tau.beta[, d, 7, setting] <- quantile(fit.2$tau.beta, probs=intervals)
    rho[, d, 7, setting] <- quantile(fit.2$rho, probs=intervals)
    nu[, d, 7, setting] <- quantile(fit.2$nu, probs=intervals)
    alpha[, d, 7, setting] <- quantile(fit.2$alpha, probs=intervals)
    
    pred <- fit.3$yp  # t5 (T=0.90)
    quant.score[, d, 8, setting] <- QuantScore(pred, probs, validate) 
    brier.score[, d, 8, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 8, setting] <- quantile(fit.3$beta[, 1], probs=intervals)
    beta.1[, d, 8, setting] <- quantile(fit.3$beta[, 2], probs=intervals)
    beta.2[, d, 8, setting] <- quantile(fit.3$beta[, 3], probs=intervals)
    tau.alpha[, d, 8, setting] <- quantile(fit.3$tau.alpha, probs=intervals)
    tau.beta[, d, 8, setting] <- quantile(fit.3$tau.beta, probs=intervals)
    rho[, d, 8, setting] <- quantile(fit.3$rho, probs=intervals)
    nu[, d, 8, setting] <- quantile(fit.3$nu, probs=intervals)
    alpha[, d, 8, setting] <- quantile(fit.3$alpha, probs=intervals)
    avgparts[, d, 3, setting] <- mean(fit.3$avgparts)
    
    pred <- fit.4$yp  # skew-t1 (T=0.90)
    quant.score[, d, 9, setting] <- QuantScore(pred, probs, validate) 
    brier.score[, d, 9, setting] <- BrierScore(pred, thresholds, validate)
    beta.0[, d, 9, setting] <- quantile(fit.4$beta[, 1], probs=intervals)
    beta.1[, d, 9, setting] <- quantile(fit.4$beta[, 2], probs=intervals)
    beta.2[, d, 9, setting] <- quantile(fit.4$beta[, 3], probs=intervals)
    tau.alpha[, d, 9, setting] <- quantile(fit.4$tau.alpha, probs=intervals)
    tau.beta[, d, 9, setting] <- quantile(fit.4$tau.beta, probs=intervals)
    rho[, d, 9, setting] <- quantile(fit.4$rho, probs=intervals)
    nu[, d, 9, setting] <- quantile(fit.4$nu, probs=intervals)
    alpha[, d, 9, setting] <- quantile(fit.4$alpha, probs=intervals)
    z.alpha[, d, 4, setting] <- quantile(fit.4$z.alpha, probs=intervals)
  }
  cat("dataset", dataset, "-b \n")	
}

save(
	quant.score, brier.score, beta.0, beta.1, beta.2, 
	tau.alpha, tau.beta, rho, nu, alpha, z.alpha, avgparts,
	probs, thresholds, 
	file="scores.RData"
)

# rm(list=ls())
# load("scores.RData")

# # dim 1: number of analysis types; dim 2: number of probs/threshs; dim 3: number of settings
# nsettings <- 9

# quant.score.mean <- array(NA, dim=c(5, length(probs), nsettings))
# brier.score.mean <- array(NA, dim=c(5, length(thresholds), nsettings))
# for(setting in 1:nsettings){
	# quant.score.mean[1,,setting] <- rowMeans(quant.score.1[,,setting])
	# quant.score.mean[2,,setting] <- rowMeans(quant.score.2[,,setting])
	# quant.score.mean[3,,setting] <- rowMeans(quant.score.3[,,setting])
	# quant.score.mean[4,,setting] <- rowMeans(quant.score.4[,,setting])
	# quant.score.mean[5,,setting] <- rowMeans(quant.score.5[,,setting])
	
	# brier.score.mean[1,,setting] <- rowMeans(brier.score.1[,,setting])
	# brier.score.mean[2,,setting] <- rowMeans(brier.score.2[,,setting])
	# brier.score.mean[3,,setting] <- rowMeans(brier.score.3[,,setting])
	# brier.score.mean[4,,setting] <- rowMeans(brier.score.4[,,setting])
	# brier.score.mean[5,,setting] <- rowMeans(brier.score.5[,,setting])
# }

# rownames(quant.score.mean) <- c("gaus", "t", "t-3", "t(.95)", "t(0.95)-3")
# colnames(quant.score.mean) <- probs
# rownames(brier.score.mean) <- c("gaus", "t", "t-3", "t(.95)", "t(0.95)-3")
# colnames(brier.score.mean) <- thresholds

# setting.title <- list()
# setting.title[[1]] <- "data: gaussian"
# setting.title[[2]] <- "data: t"
# setting.title[[3]] <- "data: t with 3 knots"
# setting.title[[4]] <- "data: 1/2 gaussian 1/2 t"
# setting.title[[5]] <- "data: 1/2 gaussian 1/2 t with 3 knots"
# setting.title[[6]] <- "data: 1/2 gaussian (range = 0.40), 1/2 t (range = 0.10)"
# setting.title[[7]] <- "data: 1/2 gaussian (range = 0.10), 1/2 t (range = 0.40)"
# setting.title[[8]] <- "data: 90% gaussian, 10% t"
# setting.title[[9]] <- "data: 95% gaussian, 5% t"

# for (setting in 1:nsettings) {
  # qsplotname <- paste("plots/all-qs-set", setting, ".pdf", sep="")
  # bsplotname <- paste("plots/all-bs-set", setting, ".pdf", sep="")
  
  # ymax <- max(quant.score.mean[, , setting])
  # ymin <- min(quant.score.mean[, , setting])
  # plot(y=quant.score.mean[1, , setting], x=probs, type='l', lty=1, ylim=c(ymin, ymax),
       # main=setting.title[[setting]], ylab="quantile scores", xaxt="n")
  
  # for (i in 2:5) {
    # lines(y=quant.score.mean[i, , setting], x=probs, lty=i)
  # }
  
  # axis(1, at=probs, labels=probs)
  # legend(x="bottomleft", legend=rownames(quant.score.mean), lty=seq(1, 5))
  # dev.print(file=qsplotname, device=pdf)

  # ymax <- max(brier.score.mean[, , setting])
  # ymin <- min(brier.score.mean[, , setting])
  # plot(y=brier.score.mean[1, , setting], x=thresholds, type='l', lty=1, ylim=c(ymin, ymax),
       # main=setting.title[[setting]], ylab="brier scores", xaxt="n")
  
  # for (i in 2:5) {
    # lines(y=brier.score.mean[i, , setting], x=thresholds, lty=i)
  # }
  
  # axis(1, at=thresholds, labels=thresholds)
  # legend(x="bottomright", legend=rownames(brier.score.mean), lty=seq(1, 5))
  # dev.print(file=bsplotname, device=pdf)

# }



# # thresholded models wrt to gaussian
# methods <- c(1, 2, 4, 5)  # removing t-3 method.
# quant.score.mean.g <- array(NA, dim=c(4, length(probs), nsettings))
# brier.score.mean.g <- array(NA, dim=c(4, length(thresholds), nsettings))

# for(i in 1:4){
	# method <- methods[i]
	# quant.score.mean.g[i, , ] <- quant.score.mean[method, , ] - quant.score.mean[1, , ]
	# brier.score.mean.g[i, , ] <- brier.score.mean[method, , ] - brier.score.mean[1, , ]
# }

# for (setting in 1:nsettings) {
  # qsplotname <- paste("plots/gaus-qs-set", setting, ".pdf", sep="")
  # bsplotname <- paste("plots/gaus-bs-set", setting, ".pdf", sep="")
  
  # ymax <- max(quant.score.mean.g[, , setting])
  # ymin <- min(quant.score.mean.g[, , setting])
  # plot(y=quant.score.mean.g[1, , setting], x=probs, type='l', lty=1, ylim=c(ymin, ymax),
       # main=setting.title[[setting]], ylab="quantile scores wrt gaussian", xaxt="n")
  
  # for (i in 2:4) {
    # lines(y=quant.score.mean.g[i, , setting], x=probs, lty=i)
  # }
  
  # axis(1, at=probs, labels=probs)
  # legend(x="topright", legend=rownames(quant.score.mean)[methods], lty=seq(1, 5))
  # dev.print(file=qsplotname, device=pdf)

  # ymax <- max(brier.score.mean.g[, , setting])
  # ymin <- min(brier.score.mean.g[, , setting])
  # plot(y=brier.score.mean.g[1, , setting], x=thresholds, type='l', lty=1, ylim=c(ymin, ymax),
       # main=setting.title[[setting]], ylab="brier scores wrt gaussian", xaxt="n")
  
  # for (i in 2:4) {
    # lines(y=brier.score.mean.g[i, , setting], x=thresholds, lty=i)
  # }
  
  # axis(1, at=thresholds, labels=thresholds)
  # legend(x="bottomleft", legend=rownames(brier.score.mean)[methods], lty=seq(1, 5))
  # dev.print(file=bsplotname, device=pdf)

# }

# # thresholded models wrt to gaussian
# methods <- c(1, 2, 4, 5)  # removing t-3 method.
# quant.score.mean.t <- array(NA, dim=c(4, length(probs), nsettings))
# brier.score.mean.t <- array(NA, dim=c(4, length(thresholds), nsettings))

# for(i in 1:4){
	# method <- methods[i]
	# quant.score.mean.t[i, , ] <- quant.score.mean[method, , ] - quant.score.mean[2, , ]
	# brier.score.mean.t[i, , ] <- brier.score.mean[method, , ] - brier.score.mean[2, , ]
# }

# for (setting in 1:nsettings) {
  # qsplotname <- paste("plots/t-qs-set", setting, ".pdf", sep="")
  # bsplotname <- paste("plots/t-bs-set", setting, ".pdf", sep="")
  
  # ymax <- max(quant.score.mean.t[, , setting])
  # ymin <- min(quant.score.mean.t[, , setting])
  # plot(y=quant.score.mean.t[1, , setting], x=probs, type='l', lty=1, ylim=c(ymin, ymax),
       # main=setting.title[[setting]], ylab="quantile scores wrt t", xaxt="n")
  
  # for (i in 2:4) {
    # lines(y=quant.score.mean.t[i, , setting], x=probs, lty=i)
  # }
  
  # axis(1, at=probs, labels=probs)
  # legend(x="topright", legend=rownames(quant.score.mean)[methods], lty=seq(1, 5))
  # dev.print(file=qsplotname, device=pdf)

  # ymax <- max(brier.score.mean.t[, , setting])
  # ymin <- min(brier.score.mean.t[, , setting])
  # plot(y=brier.score.mean.t[1, , setting], x=thresholds, type='l', lty=1, ylim=c(ymin, ymax),
       # main=setting.title[[setting]], ylab="brier scores wrt t", xaxt="n")
  
  # for (i in 2:4) {
    # lines(y=brier.score.mean.t[i, , setting], x=thresholds, lty=i)
  # }
  
  # axis(1, at=thresholds, labels=thresholds)
  # legend(x="bottomleft", legend=rownames(brier.score.mean)[methods], lty=seq(1, 5))
  # dev.print(file=bsplotname, device=pdf)

# }