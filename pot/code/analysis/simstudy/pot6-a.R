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
#  3 - skew t-1 (T = 0.90)
#  4 - skew t-3
#  5 - skew t-3 (T = 0.90)
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
source('max-stab/gen_data.R')
source('max-stab/MCMC4MaxStable.R')

setting <- 6
analysis <- "a"
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
    thresh <- quantile(y.d, probs=(0.90))
    obs <- c(rep(T, 100), rep(F, 30))
    y.o <- t(y.d[obs, ])
    x.o <- t(x[obs, , 2])
    s.o <- s[obs, ]
    knots <- s.o

    y.validate[, , d] <- t(y.d[!obs, ])
    x.p <- t(x[!obs, , 2])
    s.p <- s[!obs, ]

    cat("  start: max-stab - Set", dataset, "\n")
    tic <- proc.time()
    fit.1[[d]] <- maxstable(y=y.o, x=x.o, s=s.o, sp=s.p, xp=x.p, thresh=thresh,
                            knots=knots, iters=20000, burn=10000, update=500, thin=1)
    toc <- proc.time()
    cat("  max-stab took:", (toc - tic)[3], "\n")
    cat("  end: max-stab \n")
    cat("------------------\n")

    save(fit.1, file=outputfile)
  }
  
  rm(fit.1)
  gc()
}
