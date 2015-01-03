options(warn=2)
library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
load('us-all-setup.RData')
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

setting <- 1
method <- "gaussian"
nknots <- 1
keep.knots <- F
threshold <- 0
thresh.quant <- F
skew <- F
outputfile <- paste("OzoneFull-", setting, "US.RData", sep="")

start <- proc.time()

val <- 1

set.seed(setting*100 + val)

cat("CV", val, "started \n")
val.idx <- cv.lst[[val]]
y.o <- Y[-val.idx,]
X.o <- X[-val.idx, , ]
S.o <- S[-val.idx,]

# select every third CMAQ value
S.p <- expand.grid(x, y)
keep.these <- which((S.p[, 1] > 0) & (S.p[, 1] < 2.35) & S.p[, 2] > -1.6 & S.p[, 2] < 1.2)
CMAQ.p <- CMAQ[keep.these, ]
S.p    <- S.p[keep.these, ]

# # check the plot against the full
# zlim=c(0, 122)
# quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(CMAQ.p[, 5]), zlim=zlim)
# lines(borders/1000)

#### Thin the rows and columns by 3
# First, figure out what the x and y values are for the rows
# and columns we should keep
unique.x <- unique(S.p[, 1])
keep.x <- unique.x[seq(1, length(unique.x), by=3)]
nx <- length(keep.x)
unique.y <- unique(S.p[, 2])
keep.y <- unique.y[seq(1, length(unique(S.p[, 2])), by=3)]
ny <- length(keep.y)
keep.these <- which((S.p[, 1] %in% keep.x) & (S.p[, 2] %in% keep.y))

# Now select the subset from the original covariate and location information
CMAQ.p <- CMAQ.p[keep.these, ]
S.p    <- S.p[keep.these, ]

# # check the plot against the full
# zlim=c(0, 122)
# quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(CMAQ.p[, 5]), zlim=zlim, nx=nx, ny=ny)
# lines(borders/1000)

# reshape CMAQ to have the correct form for the mcmc function
X.p <- array(1, dim=nrow(CMAQ.p), nt, 2)
for (t in 1:nt) {
 X.p[, t, 2] <- CMAQ.p[, t]
}

tic.set <- proc.time()
fit <- mcmc(y=y.o, s=S.o, x=X.o, x.pred=X.p, s.pred=S.p,
            method=method, skew=skew, keep.knots=keep.knots,
            thresh.all=threshold, thresh.quant=thresh.quant, nknots=nknots,
            iters=30000, burn=25000, update=500, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=1,
            nu.init=0.5, alpha.init=0.5)
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

elap.time.val <- (proc.time() - start)[3]
avg.time.val <- elap.time.val / val
save(fit, file=outputfile)

