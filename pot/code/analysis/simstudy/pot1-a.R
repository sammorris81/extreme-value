#########################################################################
# A simulation study to determine if a thresholded method improves
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
# 1 - Gaussian
# 2 - skew-Gaussian (alpha = 5)
# 3 - skew-t (alpha = 5)
# 4 - skew-t w/partition (5 knots)
# 5 - 1/2 Gaussian, 1/2 t
# 6 - 1/2 Gaussian, 1/2 t w/partition (5 knots)
# 7 - 1/2 Gaussian (range = 0.40), 1/2 t (range = 0.10)
# 8 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
#
# analysis methods:
#  1 - Gaussian
#  2 - t
#  3 - t w/partition (5 knots)
#  4 - skew-Gaussian
#  5 - skew-t
#  6 - skew-t w/partition (5 knots)
#  7 - t (thresh = 0.90)
#  8 - t w/partition (5 knots, thresh = 0.90)
#  9 - skew-t (thresh=0.90)
# 10 - skew-t w/partition (5 knots, thresh=0.90)
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
analysis <- "a"
iters <- 5000; burn <- 10000; update <- 1000; thin <- 1
nsets <- 5

outputfile <- paste(setting, "-", analysis, ".RData", sep="")
set.seed(setting)
fit.1 <- fit.2 <- fit.3 <- fit.4 <- fit.5 <- vector(mode="list", length=nsets)
y.validate <- array(NA, dim=c(ntest, nt, nsets))

start <- proc.time()

for(d in 1:nsets){
  cat("start dataset", d, "\n")
  set.seed(setting * 100 + d)
  y.d <- y[, , d, setting]
  obs <- rep(c(T, F), 100)[1:ns]
  y.o <- y.d[obs, ]
  x.o <- x[obs, , ]
  s.o <- s[obs, ]

  y.validate[, , d] <- y.d[!obs, ]
  x.p <- x[!obs, , ]
  s.p <- s[!obs, ]

  cat("start: gaussian \n")
  tic <- proc.time()
  fit.1[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                     method="gaussian", skew=F, thresh=0, nknots=1,
                     iterplot=F, iters=iters, burn=burn,
                     update=update, thin=thin)
  toc <- proc.time()
  (toc - tic)[3]
  cat("end: gaussian \n")

  cat("start: t-1 \n")
  tic <- proc.time()
  fit.2[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                     method="t", skew=F, thresh=0, nknots=1,
                     iterplot=F, iters=iters, burn=burn,
                     update=update, thin=thin)
  toc <- proc.time()
  (toc - tic)[3]
  cat("end: t-1 \n")

  cat("start: t-5 \n")
  tic <- proc.time()
  fit.3[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                     method="t", skew=F, thresh=0, nknots=5,
                     iterplot=F, iters=iters, burn=burn,
                     update=update, thin=thin)
  toc <- proc.time()
  (toc - tic)[3]
  cat("end: t-5 \n")
  
  cat("start: skew gaussian \n")
  tic <- proc.time()
  fit.4[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                     method="gaussian", skew=T, thresh=0, nknots=1,
                     iterplot=F, iters=iters, burn=burn,
                     update=update, thin=thin)
  toc <- proc.time()
  (toc - tic)[3]
  cat("end: skew gaussian \n")

  cat("start: skew t-1 \n")
  tic <- proc.time()
  fit.5[[d]] <- mcmc(y=y.o, s=s.o, x=x.o, s.pred=s.p, x.pred=x.p,
                     method="t", skew=T, thresh=0, nknots=1,
                     iterplot=F, iters=iters, burn=burn,
                     update=update, thin=thin)
  toc <- proc.time()
  (toc - tic)[3]
  cat("end: skew t-1 \n")

  save(fit.1, fit.2, fit.3, fit.4, fit.5, file=outputfile)
}
