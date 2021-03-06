#########################################################################
# A simulation study to determine if a thresholded or skew methods improve
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#   1 - Gaussian
#   2 - t-1
#   3 - t-5
#   4 - skew t-1 (lambda = 3)
#   5 - skew t-5 w/partition (lambda = 3)
#   6 - max-stable with mu=1, sig=1, xi=0.2, alpha = 0.5, bw = 1
#   7 - x = setting 4, set T = q(0.80)
#       y = x,              x > T
#       y = T * exp(x - T), x <= T
#   8 - Brown-Resnick with range = 1, smooth = 0.5
#
# analysis methods:
#  1 - Gaussian
#  2 - skew t-1
#  3 - t-1 (T = 0.80)
#  4 - skew t-5
#  5 - t-5 (T = 0.80)
#  6 - Max-stable (T = 0.80)
#
#########################################################################

rm(list=ls())
load("simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- 5
ngroups <- 10
nsettings <- dim(y)[4]
nmethods <- 6
obs <- c(rep(T, 100), rep(F, 44))

setting <- 6
filename <- paste("scores", setting, ".RData", sep="")
prefix <- "results/"

# Included for easy access if we need to change score functions
source("../../R/auxfunctions.R")

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
beta.0 <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                          nsettings))
beta.1 <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                          nsettings))
beta.2 <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                          nsettings))
tau.alpha <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                             nsettings))
tau.beta  <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                             nsettings))
rho   <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                         nsettings))
nu    <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                         nsettings))
gamma <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                         nsettings))
# not all methods use skew or multiple partitions
lambda <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                          nsettings))

skew.methods <- c(2, 4)

#### Making revisions
load(filename)

quant.score.rev <- array(NA, dim=c(length(probs), (nsets * ngroups), nmethods, nsettings))
brier.score.rev <- array(NA, dim=c(length(probs), (nsets * ngroups), nmethods, nsettings))

# storage for the interval endpoints
beta.0.rev <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                              nsettings))
beta.1.rev <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                              nsettings))
beta.2.rev <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                              nsettings))
tau.alpha.rev <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                                 nsettings))
tau.beta.rev  <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                                 nsettings))
rho.rev   <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                             nsettings))
nu.rev    <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                             nsettings))
gamma.rev <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                             nsettings))
# not all methods use skew or multiple partitions
lambda.rev <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods,
                              nsettings))

quant.score.rev[, , , 1:7] <- quant.score
brier.score.rev[, , , 1:7] <- brier.score
beta.0.rev[, , , 1:7] <- beta.0
beta.1.rev[, , , 1:7] <- beta.1
beta.2.rev[, , , 1:7] <- beta.2
tau.alpha.rev[, , , 1:7] <- tau.alpha
tau.beta.rev[, , , 1:7] <- tau.beta
rho.rev[, , , 1:7] <- rho
nu.rev[, , , 1:7] <- nu
gamma.rev[, , , 1:7] <- gamma
lambda.rev[, , , 1:7] <- lambda

quant.score <- quant.score.rev
brier.score <- brier.score.rev
beta.0 <- beta.0.rev
beta.1 <- beta.1.rev
beta.2 <- beta.2.rev
tau.alpha <- tau.alpha.rev
tau.beta <- tau.beta.rev
rho <- rho.rev
nu  <- nu.rev
gamma <- gamma.rev
lambda <- lambda.rev

rm(list = c("quant.score.rev", "brier.score.rev", "beta.0.rev", "beta.1.rev",
            "beta.2.rev", "tau.alpha.rev", "tau.beta.rev", "rho.rev", "nu.rev",
            "gamma.rev", "lambda.rev"))

iters <- 20000; burn <- 10000
for (set in 1:50) {
  thresholds <- quantile(y[, , set, setting], probs=probs, na.rm=T)
  validate <- y[!obs, , set, setting]
  
  # for(setting in 1:nsettings){
  for (method in 1:nmethods) {
    if (method == 2 | method == 4) {
      dataset <- paste(prefix, setting, "-", method, "-", set, ".RData", sep="")
      load(dataset)
      
      fit <- fit.1
      if (method == 6) {
        pred <- fit$yp[10001:20000, , ]
        brier.score[, set, method, setting] <- BrierScore(pred, thresholds,
                                                          validate, trans=TRUE)
      } else {
        pred <- fit$yp
        brier.score[, set, method, setting] <- BrierScore(pred, thresholds,
                                                          validate)
        quant.score[, set, method, setting] <- QuantScore(pred, probs, validate)
      }
      
      
      if (method != 6) {
        beta.0[, set, method, setting] <- quantile(fit$beta[, 1], probs=intervals)
        beta.1[, set, method, setting] <- quantile(fit$beta[, 2], probs=intervals)
        beta.2[, set, method, setting] <- quantile(fit$beta[, 3], probs=intervals)
        tau.alpha[, set, method, setting] <- quantile(fit$tau.alpha, probs=intervals)
        tau.beta[, set, method, setting] <- quantile(fit$tau.beta, probs=intervals)
        rho[, set, method, setting] <- quantile(fit$rho, probs=intervals)
        nu[, set, method, setting] <- quantile(fit$nu, probs=intervals)
        gamma[, set, method, setting] <- quantile(fit$gamma, probs=intervals)
        if (method %in% skew.methods) {
          lambda[, set, method, setting] <- quantile(fit$lambda, probs=intervals)
        }
      }
      
      rm(fit, fit.1)
      cat("\t method", method, "\n")
    }
  }
  
  cat("dataset", set, "\n")
  
  if (set %% 10 == 0) {
    save(
      quant.score, brier.score, beta.0, beta.1, beta.2,
      tau.alpha, tau.beta, rho, nu, gamma, lambda,
      probs, thresholds,
      file = filename
    )
  }
}

