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
#   7 - x = setting 4, set T = q(0.80)
#       y = x,              x > T
#       y = T * exp(x - T), x <= T
#   8 - Brown-Resnick with range = 1, smooth = 0.5
#
# analysis methods:
#  1 - Gaussian
#  2 - skew t-1
#  3 - t-1 (T = 0.80)
#  4 - skew t-5
#  5 - t-5 (T = 0.80)
#  6 - Max-stable (T = 0.80)
#
#########################################################################

source("./package_load.R", chdir = TRUE)

# knots used in data generation
knots.x <- seq(1, 9, length=12)
knots   <- expand.grid(knots.x, knots.x)

setting <- 8
analysis <- 6
iters <- 20000; burn <- 10000; update <- 1000; thin <- 1
start <- 16; end <- 20

for (dataset in start:end) {
  outputfile <- paste(setting, "-", analysis, "-", dataset, ".RData", sep="")
  
  cat("start dataset", dataset, "\n")
  set.seed(setting * 100 + dataset)
  y.d <- y[, , dataset, setting]
  obs <- c(rep(T, 100), rep(F, 44))
  y.o <- t(y.d[obs, ])
  x.o <- x[obs, , ]  # we don't actually use this in the mcmc, we use s
  s.o <- s[obs, ]
  
  y.validate[, , d] <- y.d[!obs, ]
  x.p <- x[!obs, , ]  # we don't actually use this in the mcmc, we use s
  s.p <- s[!obs, ]
  
  thresh <- quantile(y.o, probs=0.80, na.rm=T)
  
  cat("  start: max-stable - Set", dataset, "\n")
  tic <- proc.time()
  fit.1 <- maxstable(y=y.o, x=x.o, s=s.o, sp=s.p, xp=x.p, thresh=thresh,
                     knots=knots, iters=iters, burn=burn, update=update,
                     threads=2, thin=1)
  toc <- proc.time()
  cat("  max-stable took:", (toc - tic)[3], "\n")
  cat("  end: max-stable \n")
  cat("------------------\n")
  
  save(fit.1, file=outputfile)
  rm(fit.1)
  gc()
}