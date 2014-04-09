#########################################################################
# A simulation study to determine if a thresholded method improves
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#	1 - Gaussian
#	2 - t
#	3 - t w/partition (3 knots)
#	4 - 1/2 Gaussian, 1/2 t
#	5 - 1/2 Gaussian, 1/2 t w/partition (3 knots)
#	6 - 1/2 Gaussian (range = 0.40), 1/2 t (range = 0.10)
#	7 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
#	8 - 90% Gaussian, 10% t
#	9 - 95% Gaussian, 5% t
#
# analysis methods:
#	1 - Gaussian kriging
#	2 - t
#	3 - t w/partition (3 knots)
#	4 - t (thresh = 0.95)
#	5 - t w/partition (3 knots, thresh = 0.95)
#	
#########################################################################

rm(list=ls())
load("simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsets <- dim(y)[3]
nsettings <- dim(y)[4]

# source("../../R/auxfunctions.R")	# Included for easy access if we need to change score functions

# X11()
# boxplot(log(fit.schlather$gpdre))
# lines(log(r^2))

# X11()
# lo<-apply(fit.schlather$yp,1:2,quantile,0.025)
# me<-apply(fit.schlather$yp,1:2,quantile,0.50)
# hi<-apply(fit.schlather$yp,1:2,quantile,0.975)

# plot(y.full.p,me)
# abline(0,1)
# print(mean(y.full.p>lo & y.full.p<hi))

probs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999)
thresholds <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999)


# each analysis type should have a place to score quantile scores for each set/setting
quant.score.1 <- array(NA, dim=c(length(probs), nsets, nsettings))
quant.score.2 <- array(NA, dim=c(length(probs), nsets, nsettings))
quant.score.3 <- array(NA, dim=c(length(probs), nsets, nsettings))
quant.score.4 <- array(NA, dim=c(length(probs), nsets, nsettings))
quant.score.5 <- array(NA, dim=c(length(probs), nsets, nsettings))

brier.score.1 <- array(NA, dim=c(length(thresholds), nsets, nsettings))
brier.score.2 <- array(NA, dim=c(length(thresholds), nsets, nsettings))
brier.score.3 <- array(NA, dim=c(length(thresholds), nsets, nsettings))
brier.score.4 <- array(NA, dim=c(length(thresholds), nsets, nsettings))
brier.score.5 <- array(NA, dim=c(length(thresholds), nsets, nsettings))

# complete <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
# nsettings <- length(complete)

iters <- 20000; burn <- 10000
usable <- (burn+1):iters	# only necessary if mcmc saves burnin period. ours doesn't
for(setting in 1:nsettings){
	dataset <- paste("pot-", setting, ".RData", sep="")
	load(dataset)
	for(d in 1:nsets){
		validate <- y.validate[,,d]
		quant.score.1[,d,setting] <- QuantScore(pred.1[,,,d], probs, validate) 
		quant.score.2[,d,setting] <- QuantScore(pred.2[,,,d], probs, validate)
		quant.score.3[,d,setting] <- QuantScore(pred.3[,,,d], probs, validate)
		quant.score.4[,d,setting] <- QuantScore(pred.4[,,,d], probs, validate)
		quant.score.5[,d,setting] <- QuantScore(pred.5[,,,d], probs, validate)
		
		brier.score.1[,d,setting] <- BrierScore(pred.1[,,,d], thresholds, validate)
		brier.score.2[,d,setting] <- BrierScore(pred.2[,,,d], thresholds, validate)
		brier.score.3[,d,setting] <- BrierScore(pred.3[,,,d], thresholds, validate)
		brier.score.4[,d,setting] <- BrierScore(pred.4[,,,d], thresholds, validate)
		brier.score.5[,d,setting] <- BrierScore(pred.5[,,,d], thresholds, validate)
	}
	cat("dataset", dataset)
}

save(
	quant.score.1, quant.score.2, quant.score.3, quant.score.4, quant.score.5,
	brier.score.1, brier.score.2, brier.score.3, brier.score.4, brier.score.5,
	probs, thresholds, 
	file="scores.RData"
)

rm(list=ls())
load("scores.RData")

# dim 1: number of analysis types; dim 2: number of probs/threshs; dim 3: number of settings
nsettings <- 9

quant.score.mean <- array(NA, dim=c(5, length(probs), nsettings))
brier.score.mean <- array(NA, dim=c(5, length(thresholds), nsettings))
for(setting in 1:nsettings){
	quant.score.mean[1,,setting] <- rowMeans(quant.score.1[,,setting])
	quant.score.mean[2,,setting] <- rowMeans(quant.score.2[,,setting])
	quant.score.mean[3,,setting] <- rowMeans(quant.score.3[,,setting])
	quant.score.mean[4,,setting] <- rowMeans(quant.score.4[,,setting])
	quant.score.mean[5,,setting] <- rowMeans(quant.score.5[,,setting])
	
	brier.score.mean[1,,setting] <- rowMeans(brier.score.1[,,setting])
	brier.score.mean[2,,setting] <- rowMeans(brier.score.2[,,setting])
	brier.score.mean[3,,setting] <- rowMeans(brier.score.3[,,setting])
	brier.score.mean[4,,setting] <- rowMeans(brier.score.4[,,setting])
	brier.score.mean[5,,setting] <- rowMeans(brier.score.5[,,setting])
}

rownames(quant.score.mean) <- c("gaus", "t", "t-3", "t(.95)", "t(0.95)-3")
colnames(quant.score.mean) <- probs
rownames(brier.score.mean) <- c("gaus", "t", "t-3", "t(.95)", "t(0.95)-3")
colnames(brier.score.mean) <- thresholds

setting.title <- list()
setting.title[[1]] <- "data: gaussian"
setting.title[[2]] <- "data: t"
setting.title[[3]] <- "data: t with 3 knots"
setting.title[[4]] <- "data: 1/2 gaussian 1/2 t"
setting.title[[5]] <- "data: 1/2 gaussian 1/2 t with 3 knots"
setting.title[[6]] <- "data: 1/2 gaussian (range = 0.40), 1/2 t (range = 0.10)"
setting.title[[7]] <- "data: 1/2 gaussian (range = 0.10), 1/2 t (range = 0.40)"
setting.title[[8]] <- "data: 90% gaussian, 10% t"
setting.title[[9]] <- "data: 95% gaussian, 5% t"

for (setting in 1:nsettings) {
  qsplotname <- paste("plots/all-qs-set", setting, ".pdf", sep="")
  bsplotname <- paste("plots/all-bs-set", setting, ".pdf", sep="")
  
  ymax <- max(quant.score.mean[, , setting])
  ymin <- min(quant.score.mean[, , setting])
  plot(y=quant.score.mean[1, , setting], x=probs, type='l', lty=1, ylim=c(ymin, ymax),
       main=setting.title[[setting]], ylab="quantile scores", xaxt="n")
  
  for (i in 2:5) {
    lines(y=quant.score.mean[i, , setting], x=probs, lty=i)
  }
  
  axis(1, at=probs, labels=probs)
  legend(x="bottomleft", legend=rownames(quant.score.mean), lty=seq(1, 5))
  dev.print(file=qsplotname, device=pdf)

  ymax <- max(brier.score.mean[, , setting])
  ymin <- min(brier.score.mean[, , setting])
  plot(y=brier.score.mean[1, , setting], x=thresholds, type='l', lty=1, ylim=c(ymin, ymax),
       main=setting.title[[setting]], ylab="brier scores", xaxt="n")
  
  for (i in 2:5) {
    lines(y=brier.score.mean[i, , setting], x=thresholds, lty=i)
  }
  
  axis(1, at=thresholds, labels=thresholds)
  legend(x="bottomright", legend=rownames(brier.score.mean), lty=seq(1, 5))
  dev.print(file=bsplotname, device=pdf)

}



# thresholded models wrt to gaussian
methods <- c(1, 2, 4, 5)  # removing t-3 method.
quant.score.mean.g <- array(NA, dim=c(4, length(probs), nsettings))
brier.score.mean.g <- array(NA, dim=c(4, length(thresholds), nsettings))

for(i in 1:4){
	method <- methods[i]
	quant.score.mean.g[i, , ] <- quant.score.mean[method, , ] - quant.score.mean[1, , ]
	brier.score.mean.g[i, , ] <- brier.score.mean[method, , ] - brier.score.mean[1, , ]
}

for (setting in 1:nsettings) {
  qsplotname <- paste("plots/gaus-qs-set", setting, ".pdf", sep="")
  bsplotname <- paste("plots/gaus-bs-set", setting, ".pdf", sep="")
  
  ymax <- max(quant.score.mean.g[, , setting])
  ymin <- min(quant.score.mean.g[, , setting])
  plot(y=quant.score.mean.g[1, , setting], x=probs, type='l', lty=1, ylim=c(ymin, ymax),
       main=setting.title[[setting]], ylab="quantile scores wrt gaussian", xaxt="n")
  
  for (i in 2:4) {
    lines(y=quant.score.mean.g[i, , setting], x=probs, lty=i)
  }
  
  axis(1, at=probs, labels=probs)
  legend(x="topright", legend=rownames(quant.score.mean)[methods], lty=seq(1, 5))
  dev.print(file=qsplotname, device=pdf)

  ymax <- max(brier.score.mean.g[, , setting])
  ymin <- min(brier.score.mean.g[, , setting])
  plot(y=brier.score.mean.g[1, , setting], x=thresholds, type='l', lty=1, ylim=c(ymin, ymax),
       main=setting.title[[setting]], ylab="brier scores wrt gaussian", xaxt="n")
  
  for (i in 2:4) {
    lines(y=brier.score.mean.g[i, , setting], x=thresholds, lty=i)
  }
  
  axis(1, at=thresholds, labels=thresholds)
  legend(x="bottomleft", legend=rownames(brier.score.mean)[methods], lty=seq(1, 5))
  dev.print(file=bsplotname, device=pdf)

}

# thresholded models wrt to gaussian
methods <- c(1, 2, 4, 5)  # removing t-3 method.
quant.score.mean.t <- array(NA, dim=c(4, length(probs), nsettings))
brier.score.mean.t <- array(NA, dim=c(4, length(thresholds), nsettings))

for(i in 1:4){
	method <- methods[i]
	quant.score.mean.t[i, , ] <- quant.score.mean[method, , ] - quant.score.mean[2, , ]
	brier.score.mean.t[i, , ] <- brier.score.mean[method, , ] - brier.score.mean[2, , ]
}

for (setting in 1:nsettings) {
  qsplotname <- paste("plots/t-qs-set", setting, ".pdf", sep="")
  bsplotname <- paste("plots/t-bs-set", setting, ".pdf", sep="")
  
  ymax <- max(quant.score.mean.t[, , setting])
  ymin <- min(quant.score.mean.t[, , setting])
  plot(y=quant.score.mean.t[1, , setting], x=probs, type='l', lty=1, ylim=c(ymin, ymax),
       main=setting.title[[setting]], ylab="quantile scores wrt t", xaxt="n")
  
  for (i in 2:4) {
    lines(y=quant.score.mean.t[i, , setting], x=probs, lty=i)
  }
  
  axis(1, at=probs, labels=probs)
  legend(x="topright", legend=rownames(quant.score.mean)[methods], lty=seq(1, 5))
  dev.print(file=qsplotname, device=pdf)

  ymax <- max(brier.score.mean.t[, , setting])
  ymin <- min(brier.score.mean.t[, , setting])
  plot(y=brier.score.mean.t[1, , setting], x=thresholds, type='l', lty=1, ylim=c(ymin, ymax),
       main=setting.title[[setting]], ylab="brier scores wrt t", xaxt="n")
  
  for (i in 2:4) {
    lines(y=brier.score.mean.t[i, , setting], x=thresholds, lty=i)
  }
  
  axis(1, at=thresholds, labels=thresholds)
  legend(x="bottomleft", legend=rownames(brier.score.mean)[methods], lty=seq(1, 5))
  dev.print(file=bsplotname, device=pdf)

}