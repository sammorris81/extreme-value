rm(list=ls())
load("simdata.RData")
n <- dim(y)[1]
nsets <- dim(y)[2]
nsettings <- dim(y)[3]
nthreshs <- length(thresh.probs)

# source("../../R/auxfunctions.R")	# Included for easy access if we need to change score functions

iters <- 6000; burn <- 2000
usable <- (burn + 1):iters
quant.grid <- seq(0.95, 1, 0.0001)
quant.grid <- quant.grid[-length(quant.grid)] # remove the last one for integration
grid.count <- length(quant.grid) 
grid.length <- (1 - U) / grid.count # used for integration
quant.scores <- brier.scores <- list()

for (setting in 1:nsettings) {
  dataset <- paste("thresh-", setting, ".RData", sep="")
  cat("setting", dataset, "started \n")
  load(dataset)
  quant.scores[[setting]] <- brier.scores[[setting]] <- list()
  
  for (set in 1:nsets) {
    # storage
    # quant.scores.set <- matrix(NA, nrow=grid.count, nthreshs)
    brier.scores.set <- matrix(NA, nrow=grid.count, nthreshs)
    
    for (thresh in 1:nthreshs) {
	  validate <- y[, set, setting]
	  pred.thresh <- t(preds[[set]][usable, , thresh])
	  # quant.scores.set[, thresh] <- QuantScore(preds=pred.thresh, probs=quant.grid, 
	                                           # validate=validate)
	  brier.scores.set[, thresh] <- BrierScore(preds=pred.thresh, probs=quant.grid, 
	                                           validate=validate)
    }

    # quant.scores[[setting]][[set]] <- quant.scores.set
    brier.scores[[setting]][[set]] <- brier.scores.set
    
    save(
      # quant.scores, 
      brier.scores,
      quant.grid,
      file="scores-fine.RData"
    )
    cat("    set", set, "finished \n")
  }
  
  cat("setting", dataset, "finished \n")
}


save(
    # quant.scores, 
    brier.scores,
    quant.grid,
    file="scores-fine.RData"
)

rm(list=ls())
# load("scores.RData")

# # brier scores is a list of a list of matrices
# # brier.scores[[setting]][[set]] is a matrix

# # # dim 1: number of analysis types; dim 2: number of probs/threshs; dim 3: number of settings
# U <- 0.95
# nsettings <- 3
# nsets <- 15
# grid.count <- length(quant.grid) 
# grid.length <- (1 - U) / grid.count # used for integration

# crps <- matrix(NA, nrow=3, ncol=11)

# # integrate the brier scores
# for (setting in 1:nsettings) {
  # crps.set <- matrix(NA, nrow=11, ncol=nsets)
  # for (set in 1:nsets) {
    # brier.score.set <- grid.length * brier.scores[[setting]][[set]]
    # crps.set[, set] <- colSums(brier.score.set) 
  # }
  # crps[setting, ] <- rowMeans(crps.set)
# }

# # CRPS plot
# xplot <- seq(0.85, 0.95, 0.01)
# plot(xplot, crps[1,], type="l", main="CRPS - 0.90 is true threshold")
# lines(xplot, crps[2,], lty=2)
# lines(xplot, crps[3,], lty=3)

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