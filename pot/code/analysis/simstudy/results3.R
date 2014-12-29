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
#   6 - max-stable with mu=1, sig=1, xi=0.1
#
# analysis methods:
#  1 - Gaussian
#  2 - skew t-1
#  3 - t-1 (T = 0.80)
#  4 - skew t-3
#  5 - t-3 (T = 0.80)
#  6 - max-stable
#
#########################################################################

rm(list=ls())
load("simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- 5
ngroups <- 10
done.groups <- c(7)
nsettings <- dim(y)[4]
nmethods <- 5
obs <- c(rep(T, 100), rep(F, 44))

setting <- 3
filename <- paste("scores", setting, ".RData", sep="")

source("../../R/auxfunctions.R")	# Included for easy access if we need to change score functions

# results should include
#   - coverage for all parameters
#   - quantile score plots for each data setting
#
# results do not have burnin
# fit.1[[2]] are the results for method: Gaussian on the second dataset

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)

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
z.alpha <- array(NA, dim=c(length(intervals), (nsets * ngroups), 5, nsettings))
avgparts <- array(NA, dim=c(length(intervals), (nsets * ngroups), 3, nsettings))

iters <- 20000; burn <- 10000
load(filename)
for (group in done.groups) {
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
      set.idx <- (group - 1) * 5 + d
      thresholds <- quantile(y[, , set.idx, setting], probs=probs, na.rm=T)
      validate <- y[!obs, , set.idx, setting]

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
    gc()

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
      avgparts[, set.idx, 2, setting] <- mean(fit$avgparts)
    }

    cat("dataset", dataset, "\n")
    rm(fit, fit.1)
    gc()

    # dataset <- paste(setting,"-d-",group,".RData", sep="")
    # load(dataset)
    # for(d in 1:nsets){
      # set.idx <- (group - 1) * 5 + d
      # thresholds <- quantile(y[, , set.idx, setting], probs=probs, na.rm=T)
      # validate <- y[!obs, , set.idx, setting]

      # fit <- fit.1[[d]]  # skew t-5 (T = 0.90)
      # pred <- fit$yp
      # quant.score[, set.idx, 6, setting] <- QuantScore(pred, probs, validate)
      # brier.score[, set.idx, 6, setting] <- BrierScore(pred, thresholds, validate)
      # beta.0[, set.idx, 6, setting] <- quantile(fit$beta[, 1], probs=intervals)
      # beta.1[, set.idx, 6, setting] <- quantile(fit$beta[, 2], probs=intervals)
      # beta.2[, set.idx, 6, setting] <- quantile(fit$beta[, 3], probs=intervals)
      # tau.alpha[, set.idx, 6, setting] <- quantile(fit$tau.alpha, probs=intervals)
      # tau.beta[, set.idx, 6, setting] <- quantile(fit$tau.beta, probs=intervals)
      # rho[, set.idx, 6, setting] <- quantile(fit$rho, probs=intervals)
      # nu[, set.idx, 6, setting] <- quantile(fit$nu, probs=intervals)
      # alpha[, set.idx, 6, setting] <- quantile(fit$alpha, probs=intervals)
      # z.alpha[, set.idx, 5, setting] <- quantile(fit$z.alpha, probs=intervals)
      # avgparts[, set.idx, 3, setting] <- mean(fit$avgparts)
    # }

  # # }

    # cat("dataset", dataset, "\n")
    # rm(fit, fit.1)
    # gc()

  save(
	quant.score, brier.score, beta.0, beta.1, beta.2,
	tau.alpha, tau.beta, rho, nu, alpha, z.alpha, avgparts,
	probs, thresholds,
	file=filename
  )
}
