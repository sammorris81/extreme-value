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
#  4 - skew t-5
#  5 - t-5 (T = 0.80)
#	
#########################################################################

rm(list=ls())
load("simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- 5
ngroups <- 10
nsettings <- dim(y)[4]
nmethods <- 5
obs <- rep(c(T, F), 100)[1:ns]

setting <- 1
filename <- paste("scores", setting, ".RData", sep="")

source("../../R/auxfunctions.R")	# Included for easy access if we need to change score functions

# results should include
#   - coverage for all parameters
#   - quantile score plots for each data setting
#
# results do not have burnin
# fit.1[[2]] are the results for method: Gaussian on the second dataset

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)

quant.score.all <- array(NA, dim=c(length(probs), (nsets * ngroups), nmethods, nsettings))
brier.score.all <- array(NA, dim=c(length(probs), (nsets * ngroups), nmethods, nsettings))

# storage for the interval endpoints
intervals <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
beta.0.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
beta.1.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
beta.2.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
tau.alpha.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
tau.beta.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
rho.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
nu.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
alpha.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), nmethods, nsettings))
# not all methods use skew or multiple partitions
z.alpha.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), 4, nsettings))
avgparts.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), 2, nsettings))

iters <- 20000; burn <- 10000
for (setting in 1:nsettings) {
  scores.file <- paste("scores", setting, ".RData", sep="")
  load(scores.file)
  quant.score.all[, , , setting] <- quant.score[, , , setting]
  brier.score.all[, , , setting] <- brier.score[, , , setting]
  
  beta.0.all[, , , setting] <- beta.0[, , , setting]
  beta.1.all[, , , setting] <- beta.1[, , , setting]
  beta.2.all[, , , setting] <- beta.2[, , , setting]
  tau.alpha.all[, , , setting] <- tau.alpha[, , , setting]
  rho.all[, , , setting] <- rho[, , , setting]
  nu.all[, , , setting] <- nu[, , , setting]
  alpha.all[, , , setting] <- alpha[, , , setting]
  z.alpha.all[, , , setting] <- z.alpha[, , , setting]
  avgparts.all[, , , setting] <- avgparts[, , , setting]
}  

quant.score <- quant.score.all
brier.score <- brier.score.all
beta.0 <- beta.0.all
beta.1 <- beta.1.all
beta.2 <- beta.2.all
tau.alpha <- tau.alpha.all
rho <- rho.all
nu <- nu.all
alpha <- alpha.all
z.alpha <- z.alpha.all
avgparts <- avgparts.all

rm(quant.score.all, brier.score.all, beta.0.all, beta.1.all, beta.2.all,
   tau.alpha.all, rho.all, nu.all, alpha.all, z.alpha.all, avgparts.all)

ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- 5
nsettings <- dim(y)[4]
nmethods <- 5

# get single brier scores and quantile scores for each setting x method x quantile
quant.score.mean <- apply(quant.score, c(1, 3, 4), mean, na.rm=T)
brier.score.mean <- apply(brier.score, c(1, 3, 4), mean, na.rm=T)

quant.score.mean
brier.score.mean

# find the best performing for each data setting
best.quant <- matrix(NA, nrow=11, ncol=6)
best.brier <- matrix(NA, nrow=11, ncol=6)
for (i in 1:11) {
  for (j in 1:6) {
    best.quant[i, j] <- which(quant.score.mean[i, , j] == min(quant.score.mean[i, , j]))
    best.brier[i, j] <- which(brier.score.mean[i, , j] == min(brier.score.mean[i, , j]))
  }
}

# paired t-tests
paired.results <- array(NA, dim=c(length(probs), (nmethods-1), nsettings))
compare <- c(1, 2, 4, 5)
for (i in 1:length(probs)) { for (j in 1:(nmethods-1)) { for (k in 1:nsettings) {
  if (k != 6) {
  	compare.j <- compare[j]  # want to store in the jth row of results array
  	diff <- quant.score[i, , 3, k] - quant.score[i, , compare.j, k]
  	s <- sd(diff)
  	df <- length(diff) - 1
  	t <- mean(diff) / (s / sqrt(length(diff)))
    paired.results[i, j, k] <- 2 * pt(abs(t), df=df, lower.tail=F)
  } else {
  	compare.j <- j + 1
  	diff <- quant.score[i, , 1, k] - quant.score[i, , compare.j, k]
  	s <- sd(diff)
  	df <- length(diff) - 1
  	t <- mean(diff) / (s / sqrt(length(diff)))
    paired.results[i, j, k] <- 2 * pt(abs(t), df=df, lower.tail=F)
  }
}  }  }

round(paired.results, 4)


setting.title <- c("Gaussian", "t-1", "t-5", "skew t-1 (alpha = 3)", "skew t-5 (alpha = 3)", "max-stable")
methods <- c("Gaussian", "skew-t, K = 1, T = q(0.0)", "t, K = 1, T = q(0.8)", "skew-t, K = 5, T = q(0.0)", "t, K = 5, T = q(0.8)")
bg <- c("firebrick1", "firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1")
col <- c("firebrick4", "firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4")
pch <- c(24, 22, 22, 22, 22)
lty <- c(1, 1, 1, 3, 3)

quartz(width=15, height=12)
par(mfrow=c(3, 2))
for (setting in 1:nsettings) {  
  ymax <- max(quant.score.mean[, , setting])
  ymin <- min(quant.score.mean[, , setting])
  plot(probs, quant.score.mean[, 1, setting], type='o', 
       lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
       ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]), ylab="quantile scores")
  
  for (i in 2:nmethods) {
    lines(probs, quant.score.mean[, i, setting], lty=lty[i], col=col[i])
    points(probs, quant.score.mean[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  }
  legend("topright", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg)
}

dev.print(file="plots/quantileplots.pdf", device=pdf)
dev.off()


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
