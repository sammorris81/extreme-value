load("cv.RData")
source("mcmc.R")

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

probs <- c(0.9, 0.95, 0.97, 0.99, 0.995, 0.999)
thresholds <- c(0.8, 0.85, 0.9, 0.95, 0.97, 0.99)
nsets <- 5 # Number of cv sets
nsettings <- 3 # number of different sample quantiles for thresholding
nbetas <- 4 # number of betas

quant.score.gpd <- array(NA, dim=c(length(probs), nsets, nsettings))
quant.score.gam <- array(NA, dim=c(length(probs), nsets, nsettings))
quant.score.mvn <- array(NA, dim=c(length(probs), nsets, nsettings))

brier.score.gpd <- array(NA, dim=c(length(thresholds), nsets, nsettings))
brier.score.gam <- array(NA, dim=c(length(thresholds), nsets, nsettings))
brier.score.mvn <- array(NA, dim=c(length(thresholds), nsets, nsettings))

# iters x betas x threshold x method
betas.gpd <- array(NA, dim=c(iters, nbetas, nsets, nsettings))
betas.gam <- array(NA, dim=c(iters, nbetas, nsets, nsettings))
betas.mvn <- array(NA, dim=c(iters, nbetas, nsets, nsettings)) 

usable <- (burn+1):iters
for(setting in 1:3){
	dataset <- paste("cv-", setting, ".RData", sep="")
	load(dataset)
	for(d in 1:nsets){
		val.idx <- cv.lst[[d]]
		validate <- y[val.idx,]
		quant.score.gpd[,d,setting] <- quant.score(pred.gpd[,,usable,d], probs, validate) 
		quant.score.gam[,d,setting] <- quant.score(pred.gam[,,usable,d], probs, validate)
		quant.score.mvn[,d,setting] <- quant.score(pred.mvn[,,usable,d], probs, validate)
		brier.score.gpd[,d,setting] <- brier.score(pred.gpd[,,usable,d], thresholds, validate)
		brier.score.gam[,d,setting] <- brier.score(pred.gam[,,usable,d], thresholds, validate)
		brier.score.mvn[,d,setting] <- brier.score(pred.mvn[,,usable,d], thresholds, validate)
		betas.gpd[,,d,setting] <- beta.gpd[,,d]
		betas.gam[,,d,setting] <- beta.gam[,,d]
		betas.mvn[,,d,setting] <- beta.mvn[,,d]
	}
}

savelist <- list(quant.score.gpd, quant.score.gam, quant.score.mvn, 
		     brier.score.gpd, brier.score.gam, brier.score.mvn, 
		     betas.gpd, betas.gam, betas.mvn,
		     probs, thresholds)

save(savelist, file="cv-scores.RData")
load("cv.RData")
source("mcmc.R")
load("cv-scores.RData")

quant.score.gpd <- savelist[[1]]
quant.score.gam <- savelist[[2]]
quant.score.mvn <- savelist[[3]]
brier.score.gpd <- savelist[[4]]
brier.score.gam <- savelist[[5]]
brier.score.mvn <- savelist[[6]]
betas.gpd <- savelist[[7]]
betas.gam <- savelist[[8]]
betas.mvn <- savelist[[9]]
probs <- savelist[[10]]
thresholds <- savelist[[11]]

quant.score.mean <- array(NA, dim=c(3, length(probs), 3))
brier.score.mean <- array(NA, dim=c(3, length(probs), 3))
for(setting in 1:3){
	quant.score.mean[1,,setting] <- rowMeans(quant.score.gpd[,,setting])
	quant.score.mean[2,,setting] <- rowMeans(quant.score.gam[,,setting])
	quant.score.mean[3,,setting] <- rowMeans(quant.score.mvn[,,setting])
	brier.score.mean[1,,setting] <- rowMeans(brier.score.gpd[,,setting])
	brier.score.mean[2,,setting] <- rowMeans(brier.score.gam[,,setting])
	brier.score.mean[3,,setting] <- rowMeans(brier.score.mvn[,,setting])
}

usable <- (burn+1):iters
lomehi <- c(0.025, 0.5, 0.975)
#### Setting 1: Threshold: 0.80
betas.gpd.1 <- apply(betas.gpd[usable,,,1], 2, quantile, probs=lomehi)
betas.gam.1 <- apply(betas.gam[usable,,,1], 2, quantile, probs=lomehi)
betas.mvn.1 <- apply(betas.mvn[usable,,,1], 2, quantile, probs=lomehi)

#### Setting 2: Threshold: 0.95
betas.gpd.2 <- apply(betas.gpd[usable,,,2], 2, quantile, probs=lomehi)
betas.gam.2 <- apply(betas.gam[usable,,,2], 2, quantile, probs=lomehi)
betas.mvn.2 <- apply(betas.mvn[usable,,,2], 2, quantile, probs=lomehi)

#### Setting 3: Threshold: 0.99
betas.gpd.3 <- apply(betas.gpd[usable,,,3], 2, quantile, probs=lomehi)
betas.gam.3 <- apply(betas.gam[usable,,,3], 2, quantile, probs=lomehi)
betas.mvn.3 <- apply(betas.mvn[usable,,,3], 2, quantile, probs=lomehi)