#########################################################################
# A simulation study to determine if a thresholded or skew methods improve
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#   1 - Gaussian
#   2 - t-1
#   3 - t-5
#   4 - skew t-1 (alpha = 3)
#   5 - skew t-5 w/partition (alpha = 3)
#   6 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
#
# analysis methods:
#  1 - Gaussian
#  2 - skew t-1
#  3 - skew t-1 (T = 0.90)
#  4 - skew t-5
#  5 - skew t-5 (T = 0.90)
#	
#########################################################################

rm(list=ls())
load("simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- 2
ngroups <- 10
nsettings <- dim(y)[4]
nmethods <- 5
obs <- c(rep(T, 100), rep(F, 30))

setting <- 1
filename <- paste("scores", setting, ".RData", sep="")

source("../../R/auxfunctions.R")	# Included for easy access if we need to change score functions

# results should include
#   - coverage for all parameters
#   - quantile score plots for each data setting
#
# results do not have burnin
# fit.1[[2]] are the results for method: Gaussian on the second dataset

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999)

quant.score <- array(NA, dim=c(length(probs), (nsets * ngroups), nmethods, nsettings))
brier.score <- array(NA, dim=c(length(probs), (nsets * ngroups), nmethods, nsettings))

# storage for the interval endpoints
intervals <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
beta.0 <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
beta.1 <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
beta.2 <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
tau.alpha <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
tau.beta <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
rho <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
nu <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
alpha <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
# not all methods use skew or multiple partitions
z.alpha <- array(NA, dim=c(length(intervals), (nsets * ngroups), 4, nsettings))
avgparts <- array(NA, dim=c(length(intervals), (nsets * ngroups), 2, nsettings))

iters <- 20000; burn <- 10000
group <- 1
for (group in 1) {
  # for(setting in 1:nsettings){
    dataset <- paste(setting,"-c-",group,".RData", sep="")
    load(dataset)
    
    for(d in 1:nsets){  # fit.1 is gaussian, fit.2 is t, etc.
      set.idx <- (group - 1) * 5 + d
      thresholds <- quantile(y[, , set.idx, setting], probs=probs, na.rm=T)
      validate <- y[!obs, , set.idx, setting]
    
      fit <- fit.1[[d]]  # gaussian
      pred <- fit$yp
      quant.score[, set.idx, 1, setting] <- QuantScore(pred, probs, validate) 
      brier.score[, set.idx, 1, setting] <- BrierScore(pred, thresholds, validate)
      beta.0[, set.idx, 1, setting] <- quantile(fit$beta[, 1], probs=intervals)
      beta.1[, set.idx, 1, setting] <- quantile(fit$beta[, 2], probs=intervals)
      beta.2[, set.idx, 1, setting] <- quantile(fit$beta[, 3], probs=intervals)
      tau.alpha[, set.idx, 1, setting] <- quantile(fit$tau.alpha, probs=intervals)
      tau.beta[, set.idx, 1, setting] <- quantile(fit$tau.beta, probs=intervals)
      rho[, set.idx, 1, setting] <- quantile(fit$rho, probs=intervals)
      nu[, set.idx, 1, setting] <- quantile(fit$nu, probs=intervals)
      alpha[, set.idx, 1, setting] <- quantile(fit$alpha, probs=intervals)
    
      fit <- fit.2[[d]]  # skew t-1
      pred <- fit$yp  
      quant.score[, set.idx, 2, setting] <- QuantScore(pred, probs, validate)
      brier.score[, set.idx, 2, setting] <- BrierScore(pred, thresholds, validate)
      beta.0[, set.idx, 2, setting] <- quantile(fit$beta[, 1], probs=intervals)
      beta.1[, set.idx, 2, setting] <- quantile(fit$beta[, 2], probs=intervals)
      beta.2[, set.idx, 2, setting] <- quantile(fit$beta[, 3], probs=intervals)
      tau.alpha[, set.idx, 2, setting] <- quantile(fit$tau.alpha, probs=intervals)
      tau.beta[, set.idx, 2, setting] <- quantile(fit$tau.beta, probs=intervals)
      rho[, set.idx, 2, setting] <- quantile(fit$rho, probs=intervals)
      nu[, set.idx, 2, setting] <- quantile(fit$nu, probs=intervals)
      alpha[, set.idx, 2, setting] <- quantile(fit$alpha, probs=intervals)
      z.alpha[, set.idx, 1, setting] <- quantile(fit$z.alpha, probs=intervals)
    
      fit <- fit.3[[d]]  # skew t-1 (T = 0.90)
      pred <- fit$yp  
      quant.score[, set.idx, 3, setting] <- QuantScore(pred, probs, validate)
      brier.score[, set.idx, 3, setting] <- BrierScore(pred, thresholds, validate)
      beta.0[, set.idx, 3, setting] <- quantile(fit$beta[, 1], probs=intervals)
      beta.1[, set.idx, 3, setting] <- quantile(fit$beta[, 2], probs=intervals)
      beta.2[, set.idx, 3, setting] <- quantile(fit$beta[, 3], probs=intervals)
      tau.alpha[, set.idx, 3, setting] <- quantile(fit$tau.alpha, probs=intervals)
      tau.beta[, set.idx, 3, setting] <- quantile(fit$tau.beta, probs=intervals)
      rho[, set.idx, 3, setting] <- quantile(fit$rho, probs=intervals)
      nu[, set.idx, 3, setting] <- quantile(fit$nu, probs=intervals)
      alpha[, set.idx, 3, setting] <- quantile(fit$alpha, probs=intervals)
      z.alpha[, set.idx, 2, setting] <- quantile(fit$z.alpha, probs=intervals)
    
    }
    cat("dataset", dataset, "\n")
    rm(fit, fit.1, fit.2, fit.3)
    gc()
    
    dataset <- paste(setting,"-b-",group,".RData", sep="")
    load(dataset)
    for(d in 1:nsets){ 
      fit <- fit.1[[d]]  # skew t-5
      pred <- fit$yp
      quant.score[, set.idx, 4, setting] <- QuantScore(pred, probs, validate)
      brier.score[, set.idx, 4, setting] <- BrierScore(pred, thresholds, validate)
      beta.0[, set.idx, 4, setting] <- quantile(fit$beta[, 1], probs=intervals)
      beta.1[, set.idx, 4, setting] <- quantile(fit$beta[, 2], probs=intervals)
      beta.2[, set.idx, 4, setting] <- quantile(fit$beta[, 3], probs=intervals)
      tau.alpha[, set.idx, 4, setting] <- quantile(fit$tau.alpha, probs=intervals)
      tau.beta[, set.idx, 4, setting] <- quantile(fit$tau.beta, probs=intervals)
      rho[, set.idx, 4, setting] <- quantile(fit$rho, probs=intervals)
      nu[, set.idx, 4, setting] <- quantile(fit$nu, probs=intervals)
      alpha[, set.idx, 4, setting] <- quantile(fit$alpha, probs=intervals)
      z.alpha[, set.idx, 3, setting] <- quantile(fit$z.alpha, probs=intervals)
      avgparts[, set.idx, 1, setting] <- mean(fit$avgparts)
    }
    
    cat("dataset", dataset, "\n")
    rm(fit, fit.1)
    
    dataset <- paste(setting,"-a-",group,".RData", sep="")
    load(dataset)
    for(d in 1:nsets){  
      set.idx <- (group - 1) * 5 + d
      thresholds <- quantile(y[, , set.idx, setting], probs=probs, na.rm=T)
      validate <- y[!obs, , set.idx, setting]
      
      fit <- fit.1[[d]]  # skew t-5 (T = 0.90)
      pred <- fit$yp
      quant.score[, set.idx, 5, setting] <- QuantScore(pred, probs, validate) 
      brier.score[, set.idx, 5, setting] <- BrierScore(pred, thresholds, validate)
      beta.0[, set.idx, 5, setting] <- quantile(fit$beta[, 1], probs=intervals)
      beta.1[, set.idx, 5, setting] <- quantile(fit$beta[, 2], probs=intervals)
      beta.2[, set.idx, 5, setting] <- quantile(fit$beta[, 3], probs=intervals)
      tau.alpha[, set.idx, 5, setting] <- quantile(fit$tau.alpha, probs=intervals)
      tau.beta[, set.idx, 5, setting] <- quantile(fit$tau.beta, probs=intervals)
      rho[, set.idx, 5, setting] <- quantile(fit$rho, probs=intervals)
      nu[, set.idx, 5, setting] <- quantile(fit$nu, probs=intervals)
      alpha[, set.idx, 5, setting] <- quantile(fit$alpha, probs=intervals)
      z.alpha[, set.idx, 4, setting] <- quantile(fit$z.alpha, probs=intervals)
      avgparts[, set.idx, 1, setting] <- mean(fit$avgparts)
    }
    
    cat("dataset", dataset, "\n")
    rm(fit, fit.1)
    gc()
  # }
  save(
	quant.score, brier.score, beta.0, beta.1, beta.2, 
	tau.alpha, tau.beta, rho, nu, alpha, z.alpha, avgparts,
	probs, thresholds, 
	file=filename
  )
}


# rm(list=ls())
# load("simdata.RData")
# load("scores.RData")

# ns <- dim(y)[1]
# nt <- dim(y)[2]
# nsets <- 5
# nsettings <- dim(y)[4]
# nmethods <- 5

# # get single brier scores and quantile scores for each setting x method x quantile
# quant.score.mean <- apply(quant.score, c(1, 3, 4), mean, na.rm=T)
# brier.score.mean <- apply(brier.score, c(1, 3, 4), mean, na.rm=T)

# setting.title <- c("Gaussian", "t-1", "t-5", "skew t-1 (alpha = 3)", "skew t-5 (alpha = 3)", "1/2 Gaussian (range = 0.10), 1/2 t-1 (range = 0.4)")
# methods <- c("Gaussian", "skew t-1 (T = 0.0)", "skew t-5 (T = 0.0)", "skew t-1 (T = 0.9)", "skew t-5 (T = 0.9)")
# bg <- c("firebrick1", "firebrick1", "firebrick1", "dodgerblue1", "dodgerblue1")
# col <- c("firebrick4", "firebrick4", "firebrick4", "dodgerblue4", "dodgerblue4")
# pch <- c(24, 22, 22, 22, 22)
# lty <- c(2, 1, 3, 1, 3)

# quartz(width=15, height=12)
# par(mfrow=c(3, 2))
# for (setting in 1:nsettings) {  
  # ymax <- max(quant.score.mean[, , setting])
  # ymin <- min(quant.score.mean[, , setting])
  # plot(probs, quant.score.mean[, 1, setting], type='o', 
       # lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
       # ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]), ylab="quantile scores")
  
  # for (i in 2:nmethods) {
    # lines(probs, quant.score.mean[, i, setting], lty=lty[i], col=col[i])
    # points(probs, quant.score.mean[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  # }

# }

# plot(1, 1, type="n", axes=F, main="legend", ylab="", xlab="")
# legend("center", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg, bty="n", cex=2)
# dev.print(file="plots/quantileplots.pdf", device=pdf)
# dev.off()


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
