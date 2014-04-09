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
#	6 - 1/2 Gaussian (range = 0.20), 1/2 t (range = 0.05)
#	7 - 1/2 Gaussian (range = 0.05), 1/2 t (range = 0.20)
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
source(file='../../R/mcmc.R')
source(file='../../R/auxfunctions.R')

setting <- 3
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
	r.inv.init <- r.inv.t[[setting]][, , d]
	
	y.d <- y[, , d, setting]
	obs <- rep(c(T, F), 100)[1:ns]
	y.o <- y.d[obs, ]
	x.o <- x[obs, , ]
	s.o <- s[obs, ]
	
	y.validate[, , d] <- y.d[!obs, ]
	x.p <- x[!obs, , ]
	s.p <- s[!obs, ]
	
	# iterplotname <- paste("gaus", setting, d, sep="-")
	# iterplotname <- paste(iterplotname, ".pdf", sep="")
	# cat("mvn start \n")
	# tic <- proc.time()
	# fit.1 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0, 
				  # r.model="fixed", nknots=1, iters=iters, burn=burn, 
				  # update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  # plotname=iterplotname)
	# toc <- proc.time()
	# (toc - tic)[3]
	
	# pred.1[, , , d] <- fit.1$yp[, , (burn + 1):iters]
	# beta.1[, , d]   <- fit.1$beta[(burn + 1):iters, ]
	# cat("mvn end \n")
	
	# iterplotname <- paste("t1", setting, d, sep="-")
	# iterplotname <- paste(iterplotname, ".pdf", sep="")
	# cat("t 1-knot 0-thresh start \n")
	# tic <- proc.time()
	# fit.2 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0, 
				  # r.model="gamma", nknots=1, iters=iters, burn=burn, 
				  # update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  # plotname=iterplotname)
	# toc <- proc.time()
	# (toc - tic)[3]
	
	# pred.2[, , , d] <- fit.2$yp[, , (burn + 1):iters]
	# beta.2[, , d]   <- fit.2$beta[(burn + 1):iters, ]
	# cat("t 1-knot 0-thresh end \n")
	
	# iterplotname <- paste("t3", setting, d, sep="-")
	# iterplotname <- paste(iterplotname, ".pdf", sep="")
	# cat("t 3-knots 0-thresh start \n")
	# tic <- proc.time()
	# fit.3 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0, 
				  # r.model="gamma", nknots=3, iters=iters, burn=burn, 
				  # update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  # plotname=iterplotname)
	# toc <- proc.time()
	# (toc - tic)[3]
	
	# pred.3[, , , d] <- fit.3$yp[, , (burn+1):iters]
	# beta.3[, , d]   <- fit.3$beta[(burn+1):iters, ]
	# cat("t 3-knots 0-thresh end \n")
	
	iterplotname <- paste("thresh1", setting, d, sep="-")
	iterplotname <- paste(iterplotname, ".pdf", sep="")
	cat("t 1-knot 95-thresh start \n")
	tic <- proc.time()
	fit.4 <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p, thresh=0.95, 
				  r.model="gamma", nknots=1, iters=iters, burn=burn, 
				  update=1000, thin=1, scale=T, debug=F, iterplot=T,
				  fixrho=F, fixnu=T, logrho.y.init=log(0.1), lognu.y.init=log(0.5),
				  fixbeta=F, beta.init=c(10, 0, 0),
				  fixxir=F, logxi.r.init=0, fixsigr=F, logsig.r.init=0,
				  fixr=F, r.inv.init=r.inv.init,
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
				  fixrho=T, fixnu=T, logrho.y.init=log(0.05), lognu.y.init=log(0.5),
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

y=y.o; s=s.o; x=x.o; s.pred=s.p; x.pred=x.p; thresh=0.95
r.model = "gamma" #also allow "fixed"
beta.y.m=0; beta.y.s=10; beta.init=0
logrho.y.m=-2; logrho.y.s=1; logrho.y.init=log(0.1)
lognu.y.m=-1; lognu.y.s=1; lognu.y.init=log(0.5)
alpha.y.a=1; alpha.y.b=1; alpha.y.init=0.5
logsig.r.m=0; logsig.r.s=1; logsig.r.init=0
logxi.r.m=0; logxi.r.s=0.1; logxi.r.init=0
nknots=1    
iters=20000; burn=10000; update=1000; thin=1; scale=T
debug=F; iterplot=F; plotname=NULL
fixbeta=F; fixr=F; fixrho=F; fixnu=F; fixalpha=F
fixxir=F; fixsigr=F; fixknots=F