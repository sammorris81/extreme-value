options(warn=2)
library(fields)
library(SpatialTools)
library(mvtnorm)

rm(list=ls())
load('us-all-setup.RData')
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
temporalw <- F
temporalz <- F
temporaltau <- F
outputfile <- paste("us-all-full-", setting, ".RData", sep="")

start <- proc.time()

# set.seed(4218)
# total <- 400
# NE <- S[, 1] > 0 & S[, 2] > 0
# NW <- S[, 1] < 0 & S[, 2] > 0
# SE <- S[, 1] > 0 & S[, 2] < 0
# SW <- S[, 1] < 0 & S[, 2] < 0
# nNE <- round(mean(NE) * total)
# nNW <- round(mean(NW) * total)
# nSE <- round(mean(SE) * total)
# nSW <- round(mean(SW) * total)
# observed <- c(sample(which(NE), nNE), sample(which(NW), nNW),
#                 sample(which(SE), nSE), sample(which(SW), nSW))

# set.seed(setting*100)

# y.o <- Y[observed, ]
# X.o <- X[observed, , ]
# S.o <- S[observed, ]
y.o <- Y
X.o <- X
S.o <- S


# # select every third CMAQ value
# S.p <- expand.grid(x, y)
# keep.these <- which((S.p[, 1] > -2.3) & (S.p[, 1] < 2.4) & S.p[, 2] > -1.6 & S.p[, 2] < 1.3)
# cmaq.p <- CMAQ[keep.these, ]
# S.p    <- S.p[keep.these, ]
# nx <- length(unique(S.p[, 1]))
# ny <- length(unique(S.p[, 2]))

# # # check the plot against the full
# # zlim=c(0, 122)
# # quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(CMAQ.p[, 5]), zlim=zlim, nx=nx, ny=ny)
# # lines(borders/1000)

# #### Thin the rows and columns by 3
# # First, figure out what the x and y values are for the rows
# # and columns we should keep
# unique.x <- unique(S.p[, 1])
# keep.x <- unique.x[seq(1, length(unique.x), by=10)]
# nx <- length(keep.x)
# unique.y <- unique(S.p[, 2])
# keep.y <- unique.y[seq(1, length(unique(S.p[, 2])), by=10)]
# ny <- length(keep.y)
# keep.these <- which((S.p[, 1] %in% keep.x) & (S.p[, 2] %in% keep.y))

# # Now select the subset from the original covariate and location information
# cmaq.p <- cmaq.p[keep.these, ]
# S.p    <- S.p[keep.these, ]

# # # check the plot against the full
# # dev.new()
# # zlim=c(0, 122)
# # quilt.plot(x=S.p[, 1], y=S.p[, 2], matrix(CMAQ.p[, 5]), zlim=zlim, nx=nx, ny=ny)
# # lines(borders/1000)

# # center and scale CMAQ data
# cmaq.p <- (cmaq.p - mean(cmaq.p)) / sd(cmaq.p)

# # reshape CMAQ to have the correct form for the mcmc function
# X.p <- array(1, dim=c(nrow(cmaq.p), nt, 2))
# for (t in 1:nt) {
#   X.p[, t, 2] <- cmaq.p[, t]
# }

tic.set <- proc.time()
fit <- mcmc(y=y.o, s=S.o, x=X.o, # x.pred=X.p, s.pred=S.p,
            method=method, skew=skew, keep.knots=keep.knots,
            thresh.all=threshold, thresh.quant=thresh.quant, nknots=nknots,
            iters=30000, burn=25000, update=500, iterplot=F,
            beta.init=beta.init, tau.init=tau.init, rho.init=1,
            nu.init=0.5, gamma.init=0.5)
toc.set <- proc.time()
time.set <- (toc.set - tic.set)[3]

elap.time.val <- (proc.time() - start)[3]
avg.time.val <- elap.time.val
save(fit, file=outputfile)
