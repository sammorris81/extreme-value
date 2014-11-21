load("cv-setup-us.RData")
source("../../../R/auxfunctions.R")

load("cv5-10US.RData")

par(mfrow=c(5, 5))
days1 <- c(1, 4, 7, 10, 13)
days2 <- c(16, 19, 22, 25, 28)

# plot tau
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

# plot z
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }
for(i in 6:10){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }


load("cv5-6US.RData")
par(mfrow=c(5, 5))
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$tau[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$tau[, i, days2[j]], type="l", xlab=xlab, main=main)
} }

for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days1[j]))
  plot(fit[[1]]$z[, i, days1[j]], type="l", xlab=xlab, main=main)
} }
for(i in 1:5){ for (j in 1:5) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", days2[j]))
  plot(fit[[1]]$z[, i, days2[j]], type="l", xlab=xlab, main=main)
} }


load("cv5-2US.RData")
par(mfrow=c(3, 5))
days1 <- seq(1, 30, by=2)
for (j in days1) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", j))
  plot(fit[[1]]$tau[, j], type="l", xlab=xlab, main=main)
}
for (j in days1) {
  xlab <- print(paste("knot", i))
  main <- print(paste("day", j))
  plot(fit[[1]]$z[, j], type="l", xlab=xlab, main=main)
}


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

load("cv-setup-us.RData")
source("../../../R/auxfunctions.R")

probs <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
thresholds <- quantile(Y, probs=probs, na.rm=T)
nsets <- 2 # Number of cv sets
nbetas <- 2 # number of betas

quant.score <- array(NA, dim=c(length(probs), nsets, 50))
brier.score <- array(NA, dim=c(length(thresholds), nsets, 50))

beta.0 <- array(NA, dim=c(5000, nsets, 50))
beta.1 <- array(NA, dim=c(5000, nsets, 50))

phi.z <- array(NA, dim=c(5000, nsets, 24))
phi.w <- array(NA, dim=c(5000, nsets, 24))
phi.tau <- array(NA, dim=c(5000, nsets, 24))

done <- c(1:41, 43:47, 49)

for (i in 1:50) {
  file <- paste("cv5-", i, "US.RData", sep="")
  cat("start file", file, "\n")
  if (i %in% done) {
  load(file)
  for (d in 1:nsets) {
    fit.d <- fit[[d]]
    val.idx <- cv.lst[[d]]
    validate <- Y[val.idx, ]
    pred.d <- fit.d$yp[, , ]
    if (i == 26) {
      validate <- t(validate)
      pred.d <- fit.d$yp[(25001:30000), , ]
    }
    quant.score[, d, i] <- QuantScore(pred.d, probs, validate)
    brier.score[, d, i] <- BrierScore(pred.d, thresholds, validate)
    if (i != 26) {
      beta.0[, d, i] <- fit.d$beta[, 1]
      beta.1[, d, i] <- fit.d$beta[, 2]
    }
    if (i >= 27) {
      if (i <= 42) {
        phi.z[, d, (i-26)] <- fit.d$phi.z
      }
      phi.w[, d, (i-26)] <- fit.d$phi.w
      phi.tau[, d, (i-26)] <- fit.d$phi.tau
    }
  }
  }
  cat("finish file", file, "\n")
}

savelist <- list(quant.score, brier.score,
                 beta.0, beta.1,
                 probs, thresholds)

save(savelist, file="cv-scores-us.RData")

rm(list=ls())
load("cv-setup-us.RData")
source("../../../R/auxfunctions.R")
load("cv-scores-us.RData")

quant.score <- savelist[[1]]
brier.score <- savelist[[2]]
beta.0 <- savelist[[3]]
beta.1 <- savelist[[4]]
probs <- savelist[[5]]
thresholds <- savelist[[6]]

quant.score.mean <- matrix(NA, 50, length(probs))
brier.score.mean <- matrix(NA, 50, length(thresholds))

quant.score.se <- matrix(NA, 50, length(probs))
brier.score.se <- matrix(NA, 50, length(thresholds))

done <- c(1:9, 13:41, 43, 44)
for (i in 1:50) {
  if (i %in% done) {
    quant.score.mean[i, ] <- apply(quant.score[, , i], 1, mean)
    quant.score.se[i, ] <- apply(quant.score[, , i], 1, sd) / sqrt(2)
    brier.score.mean[i, ] <- apply(brier.score[, , i], 1, mean)
    brier.score.se[i, ] <- apply(brier.score[, , i], 1, sd) / sqrt(2)
  }
}

quant.score.mean[c(1:10, 13, 14, 17:19, 21, 23, 25, 26), ]
for (i in 1:12) {
 print(which(quant.score.mean[c(1:26), i] == min(quant.score.mean[c(1:26), i])))
 print(which(brier.score.mean[c(1:26), i] == min(brier.score.mean[c(1:26), i])))
}

library(fields)
xplot <- c(1:4)
yplot <- c(1:4)
z <- matrix(quant.score.mean[c(2:17),6], nrow=4, ncol=4, byrow=F)
zlim=range(z)
image.plot(x=c(1:4), y=c(1:4), z=z, axes=F, xlab="", ylab="", main="Quantile score for q(0.95)", zlim=zlim)
text(c(row(z)), c(col(z)), "")  # not sure why we need this, but without it the labels don't turn on...
axis(side=1, at=c(1:4), labels=c("T=0", "T=50", "T=75", "T=90"), line=0.3, lwd=0)
axis(side=2, at=c(1:4), labels=c("K=1", "K=5", "K=10", "K=15"), line=0.3, las=2, lwd=0)
abline(h=c(1.5, 2.5, 3.5))
abline(v=c(1.5, 2.5, 3.5))

include.qscore <- c(1, 10, 18)
pch <- c(1, 2, 5)
plot(probs, quant.score.mean[1, ], lty=1, type="b", ylim=c(min(quant.score.mean[c(1, 10, 18), ]), max(quant.score.mean[c(1, 10, 18)])), main="Quantile Scores", xlab="quantile", ylab="score")
for (i in 2:3) {
  lines(probs, quant.score.mean[include.qscore[i], ], lty=1, type="b", pch=pch[i])
}
legend("topright", pch=c(21, 24, 23), lty=1, legend=c("Gaussian", "skew-t, K=10, T=q(0)", "Max-stable"), pt.bg="white")

round(quant.score.mean[c(2,6,10,14),c(6, 9:12)], 4)
round(quant.score.se[c(10:13,15:19),c(6, 9:12)], 4)
round(brier.score.mean[c(10:19),c(6, 9:12)]*1000, 3)
round(brier.score.se[c(10:19),c(6, 9:12)]*1000, 3)

round(quant.score.mean, 4)
round(brier.score.mean*1000, 4)

# par(mfrow=c(2, 2), oma=c(0, 0, 2, 0))

# plot(probs, quant.score.mean[1, ], lty=1, type="b", ylim=c(min(quant.score.mean), max(quant.score.mean)), main="With CMAQ", xlab="quantile", ylab="score")
# for (i in 2:9) {
  # lines(probs, quant.score.mean[i, ], lty=i)
  # points(probs, quant.score.mean[i, ], pch=i)
# }
# title(main="Quantile Scores - Ozone: NC, SC, GA", outer=T)

# plot(probs, quant.score.mean[10, ], lty=1, type="b", ylim=c(min(quant.score.mean), max(quant.score.mean)), main="Without CMAQ", xlab="quantile", ylab="score")
# for (i in 11:18) {
  # lines(probs, quant.score.mean[i, ], lty=(i-9))
  # points(probs, quant.score.mean[i, ], pch=(i-9))
# }

# plot(probs, (quant.score.mean[10, ] - quant.score.mean[1, ]), lty=1, type="b", ylim=c(0, 10), main="no CMAQ - CMAQ", xlab="quantile", ylab="score")
# for (i in 2:9) {
  # lines(probs, (quant.score.mean[(i+9), ] - quant.score.mean[i, ]), lty=i)
  # points(probs, (quant.score.mean[(i+9), ] - quant.score.mean[i, ]), pch=i)
# }

# plot(probs, quant.score.mean[1, ], type="n", axes=F)

# legend("center", lty=1:9, pch=1:9, legend=c("Gaussian", "t-1 (T=0.0)", "t-1 (T=0.9)", "t-5 (T=0.0)", "t-5 (T=0.9)", "skew-t1", "skew-t1 (T=0.9)", "skew-t5", "skew-t5 (T=0.9)"))