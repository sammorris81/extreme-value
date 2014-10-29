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
#  4 - skew t-3
#  5 - t-3 (T = 0.80)
#  6 - max-stable
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

setting <- 3
analysis <- "b"
iters <- 20000; burn <- 10000; update <- 1000; thin <- 1
nsets <- 5

for (g in 1:10) {
  fit.1 <- vector(mode="list", length=nsets)
  y.validate <- array(NA, dim=c(ntest, nt, nsets))
  outputfile <- paste(setting, "-", analysis, "-", g, ".RData", sep="")

  start <- proc.time()
  for (d in 1:nsets) {
    dataset <- (g-1) * 5 + d
    cat("start dataset", dataset, "\n")
    set.seed(setting * 100 + dataset)
    y.d <- y[, , dataset, setting]
    obs <- c(rep(T, 100), rep(F, 44))
    y.o <- y.d[obs, ]
    x.o <- x[obs, , ]
    s.o <- s[obs, ]

    y.validate[, , d] <- y.d[!obs, ]
    x.p <- x[!obs, , ]
    s.p <- s[!obs, ]

    cat("  start: skew t-5 - Set", dataset, "\n")
    tic <- proc.time()
    fit.1[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                       method="t", skew=T, thresh.all=0, thresh.quant=T,
                       nknots=5, iterplot=F, iters=iters, burn=burn,
                       update=update, thin=thin)
    toc <- proc.time()
    cat("  skew t-5 took:", (toc - tic)[3], "\n")
    cat("  end: skew t-5 \n")
    cat("------------------\n")

    save(fit.1, file=outputfile)
  }
  
  rm(fit.1)
  gc()
}
