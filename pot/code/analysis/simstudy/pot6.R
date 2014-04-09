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

library(fields)
library(geoR)
library(mvtnorm)
library(evd)

#### Load simdata
rm(list = ls())
load(file='./simdata.RData')

setting <- 6
iters <- 20000; burn <- 10000

pred.1 <- pred.2 <- pred.3 <- pred.4 <- pred.5 <- array(NA, dim=c(ntest, nt, iters-burn, nsets))
params <- array(NA, dim=c(iters, 6, nsets))
beta.1 <- beta.2 <- beta.3 <- beta.4 <- beta.5 <- array(NA, dim=c(iters-burn, 3, nsets))
y.validate <- array(NA, dim=c(ntest, nt, nsets))

outputfile <- paste("pot-", setting, ".RData", sep="")
set.seed(setting)

start <- proc.time()

for(d in 1:nsets){
	
	set.seed(setting * 100 + d)
	
	y.d <- y[, , d, setting]
	obs <- rep(c(T, F), 100)[1:ns]
	y.o <- y.d[obs, ]
	x.o <- x[obs, , ]
	s.o <- s[obs, ]
	
	y.validate[, , d] <- y.d[!obs, ]
	x.p <- x[!obs, , ]
	s.p <- s[!obs, ]
	
	iterplotname <- paste("gaus", setting, d, sep="-")
	iterplotname <- paste(iterplotname, ".pdf", sep="")
	cat("mvn start \n")
	tic <- proc.time()
	fit.1 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0, 
				  r.model="fixed", nknots=1, iters=iters, burn=burn, 
				  update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  plotname=iterplotname)
	toc <- proc.time()
	(toc - tic)[3]
	
	pred.1[, , , d] <- fit.1$yp[, , (burn + 1):iters]
	beta.1[, , d]   <- fit.1$beta[(burn + 1):iters, ]
	cat("mvn end \n")
	
	iterplotname <- paste("t1", setting, d, sep="-")
	iterplotname <- paste(iterplotname, ".pdf", sep="")
	cat("t 1-knot 0-thresh start \n")
	tic <- proc.time()
	fit.2 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0, 
				  r.model="gamma", nknots=1, iters=iters, burn=burn, 
				  update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  plotname=iterplotname)
	toc <- proc.time()
	(toc - tic)[3]
	
	pred.2[, , , d] <- fit.2$yp[, , (burn + 1):iters]
	beta.2[, , d]   <- fit.2$beta[(burn + 1):iters, ]
	cat("t 1-knot 0-thresh end \n")
	
	iterplotname <- paste("t3", setting, d, sep="-")
	iterplotname <- paste(iterplotname, ".pdf", sep="")
	cat("t 3-knots 0-thresh start \n")
	tic <- proc.time()
	fit.3 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0, 
				  r.model="gamma", nknots=3, iters=iters, burn=burn, 
				  update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  plotname=iterplotname)
	toc <- proc.time()
	(toc - tic)[3]
	
	pred.3[, , , d] <- fit.3$yp[, , (burn+1):iters]
	beta.3[, , d]   <- fit.3$beta[(burn+1):iters, ]
	cat("t 3-knots 0-thresh end \n")
	
	iterplotname <- paste("thresh1", setting, d, sep="-")
	iterplotname <- paste(iterplotname, ".pdf", sep="")
	cat("t 1-knot 95-thresh start \n")
	tic <- proc.time()
	fit.4 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0.95, 
				  r.model="gamma", nknots=1, iters=iters, burn=burn, 
				  update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  plotname=iterplotname)
	toc <- proc.time()
	(toc - tic)[3]
	
	pred.4[, , ,d] <- fit.4$yp[, , (burn + 1):iters]
	beta.4[, , d]  <- fit.4$beta[(burn+1):iters, ]
	params[, , d]  <- fit.4$params
	cat("t 1-knot 95-thresh end \n")
	
	iterplotname <- paste("thresh3", setting, d, sep="-")
	iterplotname <- paste(iterplotname, ".pdf", sep="")
	cat("t 3-knots 95-thresh start \n")
	tic <- proc.time()
	fit.5 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0.95, 
				  r.model="gamma", nknots=3, iters=iters, burn=burn, 
				  update=1000, thin=1, scale=T, debug=F, iterplot=T, 
				  plotname=iterplotname)
	toc <- proc.time()
	(toc - tic)[3]
	
	pred.5[, , , d] <- fit.5$yp[, , (burn + 1):iters]
	beta.5[, , d]   <- fit.5$beta[(burn + 1):iters, ]
	cat("t 3-knots 95-thresh end \n")	
	
	elap.time.d <- (proc.time() - start)[3]
	avg.time.d  <- elap.time.d / d
	cat("Dataset", d, "finished ", round(avg.time.d,2), " per dataset \n")
	save(pred.1, pred.2, pred.3, pred.4, pred.5,
		 beta.1, beta.2, beta.3, beta.4, beta.5,
		 params, y.validate, iters, burn, file=outputfile)
}
