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
#   6 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
#
# analysis methods:
#  1 - Gaussian
#  2 - skew t-1
#  3 - skew t-1 (T = 0.90)
#  4 - skew t-5
#  5 - skew t-5 (T = 0.90)
#	
#########################################################################

library(fields)
library(geoR)
options(warn=2)

#### Load simdata
rm(list = ls())
load(file='./simdata.RData')
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

setting <- 1
analysis <- "c"
iters <- 20000; burn <- 10000; update <- 1000; thin <- 1
nsets <- 5

for (g in 1:10) {
  fit.1 <- fit.2 <- fit.3 <- fit.4 <- fit.5 <- fit.6 <- vector(mode="list", length=nsets)
  y.validate <- array(NA, dim=c(ntest, nt, nsets))
  outputfile <- paste(setting, "-", analysis, "-", g, ".RData", sep="")

  start <- proc.time()
  for (d in 1:nsets) {
  	dataset <- (g-1) * 5 + d
    cat("start dataset", dataset, "\n")
    set.seed(setting * 100 + dataset)
    y.d <- y[, , dataset, setting]
    obs <- rep(c(T, F), 100)[1:ns]
    y.o <- y.d[obs, ]
    x.o <- x[obs, , ]
    s.o <- s[obs, ]

    y.validate[, , d] <- y.d[!obs, ]
    x.p <- x[!obs, , ]
    s.p <- s[!obs, ]

    cat("start: gaussian\n")
    tic <- proc.time()
    fit.1[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                       method="gaussian", skew=F, thresh=0, nknots=1,
                       iterplot=F, iters=iters, burn=burn,
                       update=update, thin=thin)
    toc <- proc.time()
    cat("  gaussian took:", (toc - tic)[3], "\n")
    cat("  end: gaussian \n")
    cat("------------------\n")

    cat("  start: skew t-1 \n")
    tic <- proc.time()
    fit.2[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                       method="t", skew=T, thresh=0, nknots=1,
                       iterplot=F, iters=iters, burn=burn,
                       update=update, thin=thin)
    toc <- proc.time()
    cat("  skew t-1 took:", (toc - tic)[3], "\n")
    cat("  end: skew t-1 \n")
    cat("------------------\n")

    cat("start: skew t-1 (T=0.90) \n")
    tic <- proc.time()
    fit.3[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                       method="t", skew=T, thresh=0.90, nknots=1,
                       iterplot=F, iters=iters, burn=burn,
                       update=update, thin=thin)
    toc <- proc.time()
    cat("  skew t-1 (T=0.90) took:", (toc - tic)[3], "\n")
    cat("  end: skew t-1 (T=0.90) \n")
    cat("------------------\n")

    cat("  start: skew t-5 \n")
    tic <- proc.time()
    fit.4[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                       method="t", skew=T, thresh=0, nknots=5,
                       iterplot=F, iters=iters, burn=burn,
                       update=update, thin=thin)
    toc <- proc.time()
    cat("  skew t-5 took:", (toc - tic)[3], "\n")
    cat("  end: skew t-5 \n")
    cat("------------------\n")

    cat("start: skew t-5 (T=0.90) \n")
    tic <- proc.time()
    fit.5[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                       method="t", skew=T, thresh=0.90, nknots=5,
                       iterplot=F, iters=iters, burn=burn,
                       update=update, thin=thin)
    toc <- proc.time()
    cat("  skew t-5 (T=0.90) took:", (toc - tic)[3], "\n")
    cat("  end: skew t-5 (T=0.90) \n")
    cat("------------------\n")

    save(fit.1, fit.2, fit.3, fit.4, fit.5, file=outputfile)
  }
  
  rm(fit.1, fit.2, fit.3, fit.4, fit.5)
  gc()
}
