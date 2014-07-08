#########################################################################
# A simulation study to determine if a thresholded or skew methods improve
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#	1 - Gaussian
#   2 - skew-Gaussian (alpha = 5)
#	3 - skew-t (alpha = 5)
#	4 - skew-t w/partition (5 knots)
#	5 - 1/2 Gaussian, 1/2 t
#	6 - 1/2 Gaussian, 1/2 t w/partition (5 knots)
#	7 - 1/2 Gaussian (range = 0.40), 1/2 t (range = 0.10)
#	8 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
#
# analysis methods:
#	1 - Gaussian
#   2 - skew-Gaussian
#   3 - skew-t
#   4 - skew-t w/partition (5 knots)
#	5 - t
#	6 - t w/partition (5 knots)
#	7 - t (thresh = 0.90)
#	8 - t w/partition (5 knots, thresh = 0.90)
#	
#########################################################################

# clean out previous variables
rm(list=ls())

# libraries
library(fields)
library(geoR)

# necessary functions
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

# data settings
beta.t <- c(10, 2, -3)
nu.t <- 0.5
alpha.t <- 0.9
mixprob.t <- c(0, 1, 1, 1, 0.5, 0.5, 0.5, 0.5)  # 0: Gaussian, 1: t
nknots.t <- c(1, 1, 1, 5, 1, 5, 1, 1)
gau.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.40, 0.10)
t.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.40)
z.alpha.t <- c(0, 3, 3, 3, 0, 0, 0, 0)
tau.alpha.t <- 2
tau.beta.t  <- 8


# covariate data
s         <- cbind(runif(60), runif(60))
ns        <- nrow(s)  
nt        <- 50
nsets     <- 20
nsettings <- length(mixprob.t) 
ntest     <- floor(ns / 2)

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

# Storage for datasets
y <- array(NA, dim=c(ns, nt, nsets, nsettings))
tau.t <- z.t <- knots.t <- vector("list", length=nsettings)

for (setting in 1:nsettings) {
  nknots <- nknots.t[setting]
  tau.t.setting <- array(NA, dim=c(nknots, nt, nsets))
  z.t.setting <- array(NA, dim=c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim=c(nknots, 2, nt, nsets))
  for (set in 1:nsets) {
    set.seed(setting*100 + set)
    data <- rpotspat(nt=nt, x=x, s=s, beta=beta.t, alpha=alpha.t, nu=nu.t,
                     gau.rho=gau.rho.t[setting], t.rho=t.rho.t[setting],
                     mixprob=mixprob.t[setting], z.alpha=z.alpha.t[setting],
                     tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                     nknots=nknots.t[setting])
    
    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }
    
  tau.t[[setting]]   <- tau.t.setting
  z.t[[setting]]     <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
}

# remove simstudy "truth" settings to help diagnose errors.
save(y, tau.t, z.t, knots.t, ns, nt, s, nsets, 
     x, ntest,  # covariate data that should be the same for all datasets
     file='simdata.RData')
rm(list=ls())
