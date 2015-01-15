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
done.groups <- c(1:6)
done.sets <- rep(NA, length(done.groups) * nsets)
for (i in 1:length(done.groups)){
  idx <- (i - 1) * 5 + seq(1:5)
  done.sets[idx] <- (done.groups[i] - 1) * 5 + seq(1:5)
}
nsettings <- 6
nmethods <- 5
obs <- c(rep(T, 100), rep(F, 44))

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
z.alpha.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), 5, nsettings))
avgparts.all <- array(NA, dim=c(length(intervals), (nsets * ngroups), 3, nsettings))

iters <- 20000; burn <- 10000
for (setting in 1:nsettings) {
  scores.file <- paste("scores", setting, ".RData", sep="")
  load(scores.file)
  quant.score.all[, done.sets, , setting] <- quant.score[, done.sets, , setting]
  brier.score.all[, done.sets, , setting] <- brier.score[, done.sets, , setting]

  beta.0.all[, done.sets, , setting] <- beta.0[, done.sets, , setting]
  beta.1.all[, done.sets, , setting] <- beta.1[, done.sets, , setting]
  beta.2.all[, done.sets, , setting] <- beta.2[, done.sets, , setting]
  tau.alpha.all[, done.sets, , setting] <- tau.alpha[, done.sets, , setting]
  rho.all[, done.sets, , setting] <- rho[, done.sets, , setting]
  nu.all[, done.sets, , setting] <- nu[, done.sets, , setting]
  alpha.all[, done.sets, , setting] <- alpha[, done.sets, , setting]
  z.alpha.all[, done.sets, , setting] <- z.alpha[, done.sets, , setting]
  avgparts.all[, done.sets, , setting] <- avgparts[, done.sets, , setting]
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
nsettings <- 6
nmethods <- 5

# get single brier scores and quantile scores for each setting x method x quantile
quant.score.mean <- apply(quant.score, c(1, 3, 4), mean, na.rm=T)
brier.score.mean <- apply(brier.score, c(1, 3, 4), mean, na.rm=T)

quant.score.mean
brier.score.mean

# find the best performing for each data setting
best.quant <- matrix(NA, nrow=11, ncol=nsettings)
best.brier <- matrix(NA, nrow=11, ncol=nsettings)
for (i in 1:11) {
  for (j in 1:nsettings) {
    best.quant[i, j] <- which(quant.score.mean[i, , j] == min(quant.score.mean[i, , j], na.rm=T))
    best.brier[i, j] <- which(brier.score.mean[i, , j] == min(brier.score.mean[i, , j], na.rm=T))
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
wilcox.results.gau <- array(NA, dim=c(length(probs), (nmethods-1), nsettings))
for (i in 1:length(probs)) { for (k in 1:nsettings) {
  ref.j <- 1
  for (j in 2:nmethods) {
    wilcox.results.gau[i, (j-1), k] <- wilcox.test(brier.score[i, , ref.j, k], brier.score[i, , j, k],
                                               paired=T, alternative="greater")$p.value
  }
}}

wilcox.results.5 <- matrix(NA, nrow=length(probs), ncol=nsettings)
for (i in 1:length(probs)) { for (k in 1:nsettings) {
  wilcox.results.5[i, k] <- wilcox.test(brier.score[i, , 2, k], brier.score[i, , 4, k],
                                               paired=T)$p.value
}}

# get single brier scores and quantile scores for each setting x method x quantile
quant.score.med <- apply(quant.score, c(1, 3, 4), median, na.rm=T)
brier.score.med <- apply(brier.score, c(1, 3, 4), median, na.rm=T)

bs.med.ref.gau <- array(NA, dim=c(11, 4, nsettings))
bs.mean.ref.gau <- array(NA, dim=c(11, 4, nsettings))
for (j in 1:4) {
  bs.med.ref.gau[, j, ] <- brier.score.med[, (j + 1), ] / brier.score.med[, 1, ]
  bs.mean.ref.gau[, j, ] <- brier.score.mean[, (j + 1), ] / brier.score.mean[, 1, ]
}

qs.med.ref.gau <- array(NA, dim=c(11, 4, nsettings))
qs.mean.ref.gau <- array(NA, dim=c(11, 4, nsettings))
for (j in 1:4) {
  qs.med.ref.gau[, j, ] <- quant.score.med[, (j + 1), ] / quant.score.med[, 1, ]
  qs.mean.ref.gau[, j, ] <- quant.score.mean[, (j + 1), ] / quant.score.mean[, 1, ]
}


setting.title <- c("Gaussian", "Symmetric-t (K = 1)", "Symmetric-t (K = 5)", "Skew-t (K = 1, alpha = 3)", "Skew-t (K = 5, alpha = 3)", "Max-stable")
methods <- c("Skew-t, K = 1, T = q(0.0)", "Skew-t, K = 1, T = q(0.8)", "Skew-t, K = 5, T = q(0.0)", "Skew-t, K = 5, T = q(0.8)")
bg <- c("firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1")
col <- c("firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4")
pch <- c(22, 22, 22, 22)
lty <- c(1, 1, 3, 3)

# Individual plots for presentation
setting <- 1
ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T)
ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
plot(probs, bs.mean.ref.gau[, 1, setting], type='o',
     lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
     ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]),
     ylab="Relative Brier score", xlab="Threshold quantile")

for (i in 2:(nmethods - 1)) {
  lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
  points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  abline(h=1, lty=2)
}

setting <- 2
ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T)
ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
plot(probs, bs.mean.ref.gau[, 1, setting], type='o',
     lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
     ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]),
     ylab="Relative Brier score", xlab="Threshold quantile")

for (i in 2:(nmethods - 1)) {
  lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
  points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  abline(h=1, lty=2)
}

setting <- 3
ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T)
ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
plot(probs, bs.mean.ref.gau[, 1, setting], type='o',
     lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
     ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]),
     ylab="Relative Brier score", xlab="Threshold quantile")

for (i in 2:(nmethods - 1)) {
  lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
  points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  abline(h=1, lty=2)
}

setting <- 4
ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T)
ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
plot(probs, bs.mean.ref.gau[, 1, setting], type='o',
     lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
     ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]),
     ylab="Relative Brier score", xlab="Threshold quantile")

for (i in 2:(nmethods - 1)) {
  lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
  points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  abline(h=1, lty=2)
}

setting <- 5
ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T)
ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
plot(probs, bs.mean.ref.gau[, 1, setting], type='o',
     lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
     ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]),
     ylab="Relative Brier score", xlab="Threshold quantile")

for (i in 2:(nmethods - 1)) {
  lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
  points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  abline(h=1, lty=2)
}

setting <- 6
ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T)
ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
plot(probs, bs.mean.ref.gau[, 1, setting], type='o',
     lty=lty[1], pch=pch[1], col=col[1], bg=bg[1],
     ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]),
     ylab="Relative Brier score", xlab="Threshold quantile")

for (i in 2:(nmethods - 1)) {
  lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
  points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i])
  abline(h=1, lty=2)
}


plot(1, 1, type='n', axes=F, ylab="", xlab="")
legend("center", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg)


# Panel for paper
quartz(width=15, height=12)
par(mfrow=c(3, 2), mar=c(5.1, 5.1, 4.1, 2.1))
for (setting in 1:nsettings) {
  if (setting == 6) {
    ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T) + 0.1
    ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
  } else {
    ymax <- max(bs.mean.ref.gau[, , setting], 1, na.rm=T)
    ymin <- min(bs.mean.ref.gau[, , setting], 1, na.rm=T)
  }

  plot(probs, bs.mean.ref.gau[, 1, setting], type='o',
       lty=lty[1], pch=pch[1], col=col[1], bg=bg[1], cex=1.5,
       ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]), ylab="Relative Brier score", xlab="Threshold quantile", cex.lab=2, cex.axis=2, cex.main=2)

  for (i in 2:(nmethods - 1)) {
    lines(probs, bs.mean.ref.gau[, i, setting], lty=lty[i], col=col[i])
    points(probs, bs.mean.ref.gau[, i, setting], pch=pch[i], col=col[i], bg=bg[i], cex=1.5)
    abline(h=1, lty=2)
  }
  if (setting == 6) {
  	legend("topright", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg, cex=1.75)
  }
}

dev.print(file="plots/bsplots-mean.pdf", device=pdf)
dev.off()

setting.title <- c("Gaussian", "T (K = 1)", "T (K = 5)", "Skew-t (K = 1, alpha = 3)", "Skew-t (K = 5, alpha = 3)", "Max-stable")
methods <- c("Skew-t, K = 1, T = q(0.0)", "Skew-t, K = 1, T = q(0.8)", "Skew-t, K = 5, T = q(0.0)", "Skew-t, K = 5, T = q(0.8)")
bg <- c("firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1")
col <- c("firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4")
pch <- c(22, 22, 22, 22)
lty <- c(1, 1, 3, 3)

quartz(width=15, height=12)
par(mfrow=c(3, 2), mar=c(5.1, 5.1, 4.1, 2.1))
for (setting in 1:nsettings) {
  if (setting == 6) {
    ymax <- max(qs.mean.ref.gau[-c(10,11), , setting], 1) + 0.05
  } else {
    ymax <- max(qs.mean.ref.gau[-c(10,11), , setting], 1)
  }
  ymin <- min(qs.mean.ref.gau[-c(10,11), , setting], 1)
  plot(probs[-c(10,11)], qs.mean.ref.gau[-c(10,11), 1, setting], type='o',
       lty=lty[1], pch=pch[1], col=col[1], bg=bg[1], cex=1.5,
       ylim=c(ymin, ymax), main=paste("Data:", setting.title[setting]), ylab="Relative quantile score", xlab="Threshold quantile", cex.lab=2, cex.axis=2, cex.main=2)

  for (i in 2:(nmethods - 1)) {
    lines(probs[-c(10,11)], qs.mean.ref.gau[-c(10,11), i, setting], lty=lty[i], col=col[i])
    points(probs[-c(10,11)], qs.mean.ref.gau[-c(10,11), i, setting], pch=pch[i], col=col[i], bg=bg[i], cex=1.5)
    abline(h=1, lty=2)
  }
  if (setting == 6) {
    legend("topleft", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg, cex=1.6)
  }
}

dev.print(file="plots/qsplots-mean.pdf", device=pdf)
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