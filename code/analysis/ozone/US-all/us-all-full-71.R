options(warn=2)
library(fields)
library(SpatialTools)
library(mvtnorm)

rm(list=ls())
load("../ozone_data.RData")
source('../../../R/mcmc.R', chdir=T)
source('../../../R/auxfunctions.R')

setting <- 71
method <- "t"
nknots <- 10
keep.knots <- TRUE
threshold <- 75
tau.init <- 0.05
thresh.quant <- FALSE
skew <- FALSE
temporalw <- TRUE
temporalz <- TRUE
temporaltau <- TRUE
beta.init <- 0
tau.init <- 1
outputfile <- paste("us-all-full-", setting, ".RData", sep="")

# rescale x and y coordinates to make easier to work with
x <- x / 1000
y <- y / 1000

S       <- cbind(x[s[, 1]], y[s[, 2]])  # expands the grid of x, y
excl    <- which(rowMeans(is.na(Y)) > 0.50)  # remove where we're missing 50%
index   <- index[-excl]
Y       <- Y[-excl, ]
S       <- S[-excl, ]
CMAQ.cs <- (CMAQ - mean(CMAQ)) / sd(CMAQ)  # center and scale CMAQ data
cmaq    <- CMAQ.cs[index, ]  # extract cmaq for sites

# make design matrix
nt <- ncol(Y)
X  <- array(1, dim=c(nrow(cmaq), nt, 2))
for (t in 1:nt) {
  X[, t, 2] <- cmaq[, t]
}

start <- proc.time()

y.o <- Y
X.o <- X
S.o <- S

tic.set <- proc.time()
fit <- mcmc(y=y.o, s=S.o, x=X.o, # x.pred=X.p, s.pred=S.p,
            method=method, skew=skew, keep.knots=keep.knots,
            min.s=c(-2.25, -1.60), max.s=c(2.35, 1.30),
            thresh.all=threshold, thresh.quant=thresh.quant, nknots=nknots,
            iters=30000, burn=25000, update=500, iterplot=F,
            beta.init=beta.init, tau.init=tau.init,
            gamma.init=0.5, rho.init=1, rho.upper=5, nu.init=0.5, nu.upper=10,
            temporaltau=temporaltau, temporalw=temporalw, temporalz=temporalz)
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

elap.time.val <- (proc.time() - start)[3]
avg.time.val <- elap.time.val
save(fit, file=outputfile)
