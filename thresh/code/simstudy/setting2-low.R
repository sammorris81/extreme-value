library(evd)
rm(list=ls())
load(file='simdata.RData')

setting <- 2

n <- dim(y)[1]
nsets <- dim(y)[2]
iters <- 6000; burn <- 2000
thresh.probs <- seq(0.50, 0.95, 0.05)

nthreshs <- length(thresh.probs)
outputfile <- paste("thresh-low", setting, ".RData", sep="")


# want to store all of the predicted datasets - list that is nsets long
preds <- beta.prob <- beta.scale <- beta.shape <- list()

for (set in 1:nsets) {
  cat("Dataset", set, "started \n")
  tic <- proc.time()
  y.set <- y[, set, setting]
  
  threshs <- quantile(y.set, probs=thresh.probs)
  
  preds[[set]] <- array(NA, dim=c(iters, n, nthreshs))
  beta.prob[[set]] <- beta.scale[[set]] <- beta.shape[[set]] <- array(NA, dim=c(iters, ncol(X), nthreshs))
  
  for (t in 1:nthreshs){
    fit <- Bayes_GPD(y=y.set, X.prob=X, X.sig=X, X.xi=X, 
                     Xp.prob=X, Xp.sig=X, Xp.xi=X,
                     thresh=threshs[t], iters=iters, burn=burn, update=1000,
                     iterplot=F)
  
    beta.prob[[set]][, , t] <- fit$beta.prob
    beta.scale[[set]][, , t] <- fit$beta.sig
    beta.shape[[set]][, , t] <- fit$beta.xi
    preds[[set]][, , t] <- fit$yp
  	cat("  thresh", t, "finished \n")
  	
  	save(beta.prob, beta.scale, beta.shape, preds,
  	     file=outputfile)
  }
  
  elap.time.set <- (proc.time() - tic)[3]
  avg.time.set  <- elap.time.set / set
  cat("Dataset", set, "finished ", round(avg.time.set, 2), "per dataset \n")
}
