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
# Look at coverage for non-thresholded methods that match data setting
rm(list=ls())
load("simdata.RData")
beta.t <- c(10, 0, 0)
nu.t <- 0.5
alpha.t <- 0.9
mixprob.t <- c(0, 1, 1, 1, 1, 0.5)  # 0: Gaussian, 1: t
nknots.t <- c(1, 1, 5, 1, 5, 1)
gau.rho.t <- c(1, 1, 1, 1, 1, 1)
t.rho.t <- c(1, 1, 1, 1, 1, 4)
z.alpha.t <- c(0, 0, 0, 3, 3, 0)
tau.alpha.t <- 3
tau.beta.t  <- 8

# load in the dataset
grp <- 1
filename <- paste("2-c-", grp, ".RData", sep="")
load(filename)
par(mfrow=c(3, 5))
for (i in 1:5) {
  interval <- quantile(fit.2[[i]]$tau.alpha, probs=c(0.025, 0.975))
  plot(fit.2[[i]]$tau.alpha, type="l", ylab="tau.alpha", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (i in 1:5) {
  interval <- round(quantile(fit.2[[i]]$tau.beta, probs=c(0.025, 0.975)), 1)
  plot(fit.2[[i]]$tau.beta, type="l", ylab="tau.beta", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))  
}
for (i in 1:5) {
  interval <- round(quantile(fit.2[[i]]$z.alpha, probs=c(0.025, 0.975)), 1)
  plot(fit.2[[i]]$z.alpha, type="l", ylab="z.alpha", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}

par(mfrow=c(3, 5))
for (i in 1:5) {
  interval <- round(quantile(fit.2[[i]]$rho, probs=c(0.025, 0.975)), 2)
  plot(fit.2[[i]]$rho, type="l", ylab="rho", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (i in 1:5) {
  interval <- round(quantile(fit.2[[i]]$nu, probs=c(0.025, 0.975)), 2)
  plot(fit.2[[i]]$nu, type="l", ylab="nu", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))  
}
for (i in 1:5) {
  interval <- round(quantile(fit.2[[i]]$alpha, probs=c(0.025, 0.975)), 2)
  plot(fit.2[[i]]$alpha, type="l", ylab="alpha", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}

par(mfrow=c(5, 5))
daygrp <- 1
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.2[[1]]$tau[, day], probs=c(0.025, 0.975)), 2)
  plot(fit.2[[1]]$tau[, day], type="l", ylab=round(tau.t[[2]][1, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
daygrp <- 2
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.2[[1]]$tau[, day], probs=c(0.025, 0.975)), 2)
  plot(fit.2[[1]]$tau[, day], type="l", ylab=round(tau.t[[2]][1, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}


grp <- 1
filename <- paste("3-b-", grp, ".RData", sep="")
load(filename)
par(mfrow=c(3, 5))
for (i in 1:5) {
  interval <- quantile(fit.1[[i]]$tau.alpha, probs=c(0.025, 0.975))
  plot(fit.1[[i]]$tau.alpha, type="l", ylab="tau.alpha", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (i in 1:5) {
  interval <- round(quantile(fit.1[[i]]$tau.beta, probs=c(0.025, 0.975)), 1)
  plot(fit.1[[i]]$tau.beta, type="l", ylab="tau.beta", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))  
}
for (i in 1:5) {
  interval <- round(quantile(fit.1[[i]]$z.alpha, probs=c(0.025, 0.975)), 1)
  plot(fit.1[[i]]$z.alpha, type="l", ylab="z.alpha", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}

par(mfrow=c(3, 5))
for (i in 1:5) {
  interval <- round(quantile(fit.1[[i]]$rho, probs=c(0.025, 0.975)), 2)
  plot(fit.1[[i]]$rho, type="l", ylab="rho", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (i in 1:5) {
  interval <- round(quantile(fit.1[[i]]$nu, probs=c(0.025, 0.975)), 2)
  plot(fit.1[[i]]$nu, type="l", ylab="nu", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))  
}
for (i in 1:5) {
  interval <- round(quantile(fit.1[[i]]$alpha, probs=c(0.025, 0.975)), 2)
  plot(fit.1[[i]]$alpha, type="l", ylab="alpha", xlab=print(paste("Set:", (grp - 1) * 5 +i)),
  main=print(paste(interval[1], ", ", interval[2])))
}

par(mfrow=c(5, 5))
daygrp <- 1
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 1, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 1, day], type="l", ylab=round(tau.t[[3]][1, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}

for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 2, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 2, day], type="l", ylab=round(tau.t[[3]][2, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 3, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 3, day], type="l", ylab=round(tau.t[[3]][3, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 4, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 4, day], type="l", ylab=round(tau.t[[3]][4, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 5, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 5, day], type="l", ylab=round(tau.t[[3]][5, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}

daygrp <- 2
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 1, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 1, day], type="l", ylab=round(tau.t[[3]][1, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 2, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 2, day], type="l", ylab=round(tau.t[[3]][2, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 3, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 3, day], type="l", ylab=round(tau.t[[3]][3, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 4, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 4, day], type="l", ylab=round(tau.t[[3]][4, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}
for (t in 1:25) {
  day <- (daygrp - 1) * 25 + t
  interval <- round(quantile(fit.1[[1]]$tau[, 5, day], probs=c(0.025, 0.975)), 2)
  plot(fit.1[[1]]$tau[, 5, day], type="l", ylab=round(tau.t[[3]][5, day, 1], 2), xlab=print(paste("Day:", day)),
  main=print(paste(interval[1], ", ", interval[2])))
}

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

# wilcoxon signed-rank test for Brier scores
wilcox.results <- array(NA, dim=c(length(probs), (nmethods-1), nsettings))
best.methods <- vector(mode="list")
for (i in 1:length(probs)) { for (k in 1:nsettings) {
  ref.j <- best.brier[i, k]
  cmp.j <- (1:5)[-ref.j]
  for (j in 1:(nmethods-1)) {
    wilcox.results[i, j, k] <- wilcox.test(brier.score[i, , ref.j, k], brier.score[i, , cmp.j[j], k], paired=T)$p.value
  }
}}

# get single brier scores and quantile scores for each setting x method x quantile
quant.score.med <- apply(quant.score, c(1, 3, 4), median, na.rm=T)
brier.score.med <- apply(brier.score, c(1, 3, 4), median, na.rm=T)

bs.med.ref.gau <- array(NA, dim=c(11, 4, 6))
bs.mean.ref.gau <- array(NA, dim=c(11, 4, 6))
for (j in 1:4) {
  bs.med.ref.gau[, j, ] <- brier.score.med[, (j + 1), ] / brier.score.med[, 1, ]
  bs.mean.ref.gau[, j, ] <- brier.score.mean[, (j + 1), ] / brier.score.mean[, 1, ]
}


setting.title <- c("Gaussian", "t (K = 1)", "t (K = 5)", "skew t (K = 1, alpha = 3)", "skew t (K = 5, alpha = 3)", "max-stable")
methods <- c("skew-t, K = 1, T = q(0.0)", "t, K = 1, T = q(0.8)", "skew-t, K = 5, T = q(0.0)", "t, K = 5, T = q(0.8)")
bg <- c("firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1")
col <- c("firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4")
pch <- c(22, 22, 22, 22)
lty <- c(1, 1, 3, 3)

quartz(width=15, height=12)
par(mfrow=c(3, 2))
for (setting in 1:nsettings) {  
  ymax <- max(bs.med.ref.gau[, , setting], 1)
  ymin <- min(bs.med.ref.gau[, , setting], 1)
  plot(probs, bs.med.ref.gau[, 1, setting], type='o', 
       lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
       ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]), ylab="relative brier score")
  
  for (i in 2:(nmethods - 1)) {
    lines(probs, bs.med.ref.gau[, i, setting], lty=lty[i], col=col[i])
    points(probs, bs.med.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
    abline(h=1, lty=3)
  }
  if (setting == 6) {
  	legend("bottomleft", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg, cex=1.5)
  }
}

dev.print(file="plots/bsplots-med.pdf", device=pdf)
dev.off()

setting.title <- c("Gaussian", "T (K = 1)", "T (K = 5)", "Skew-t (K = 1, alpha = 3)", "Skew-t (K = 5, alpha = 3)", "Max-stable")
methods <- c("Skew-t, K = 1, T = q(0.0)", "T, K = 1, T = q(0.8)", "Skew-t, K = 5, T = q(0.0)", "T, K = 5, T = q(0.8)")
bg <- c("firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1")
col <- c("firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4")
pch <- c(22, 22, 22, 22)
lty <- c(1, 1, 3, 3)

quartz(width=15, height=12)
par(mfrow=c(3, 2), mar=c(5.1, 5.1, 4.1, 2.1))
for (setting in 1:nsettings) {  
  ymax <- max(bs.mean.ref.gau[, , setting], 1)
  ymin <- min(bs.mean.ref.gau[, , setting], 1)
  plot(probs, bs.mean.ref.gau[, 1, setting], type='o', 
       lty=lty[1], pch=pch[1], col=col[1], bg=bg[1], cex=1.5,
       ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]), ylab="Relative brier score", xlab="Threshold quantile", cex.lab=2, cex.axis=2, cex.main=2)
  
  for (i in 2:(nmethods - 1)) {
    lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
    points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i], cex=1.7)
    abline(h=1, lty=2)
  }
  if (setting == 6) {
  	legend("bottomleft", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg, cex=1.7)
  }
}

dev.print(file="plots/bsplots-mean.pdf", device=pdf)
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
