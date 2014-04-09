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
#	6 - 1/2 Gaussian (range = 0.40), 1/2 t (range = 0.10)
#	7 - 1/2 Gaussian (range = 0.10), 1/2 t (range = 0.40)
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

# clean out previous variables
rm(list=ls())

# libraries
library(fields)
library(geoR)
library(mvtnorm)
library(evd)

# necessary functions
# setwd('~/repos-svn/dissertation/extreme-value/potv2/code/analysis/simstudy')
source('../../R/mcmc.R')
source('../../R/densities.R')
source('../../R/auxfunctions.R')

# data settings
beta.y.t <- c(10, 0, 0)
sig2.y.t <- 1
nu.y.t <- 0.5
alpha.y.t <- 0.9
mixprob.t <- c(1, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.9, 0.95)	# 1: r=1; 0: r=IG
nknots.t <- c(1, 1, 3, 1, 3, 1, 1, 1, 1)
mvn.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.40, 0.10, 0.10, 0.10)
gam.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.40, 0.10, 0.10)
xi.r.t <- 1     # first parameter for IG
sig.r.t <- 1    # second paratmeter for IG 
thresh  <- 0.90

# covariate data
s         <- cbind(runif(50), runif(50))
ns        <- nrow(s)  
nt        <- 10
nsets     <- 15
nsettings <- length(mixprob.t) 
ntest     <- floor(ns / 2)

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

# Storage for datasets
y <- array(NA, dim=c(ns, nt, nsets, nsettings))
r.inv.t <- knots.t <- vector("list", length=nsettings)

for (setting in 1:nsettings) {
	nknots <- nknots.t[setting]
	r.inv.t.setting <- array(NA, dim=c(nknots, nt, nsets))
	knots.t.setting <- array(NA, dim=c(nknots, 2, nt, nsets))
    for (set in 1:nsets) {
    set.seed(setting*100 + set)
    data <- rpotspatial(nt=nt, s=s, x=x, beta.y=beta.y.t, nu.y=nu.y.t, 
                        alpha.y=alpha.y.t, thresh=thresh, xi.r=xi.r.t, 
                        sig.r=sig.r.t, mixprob=mixprob.t[setting],
                        mvn.rho=mvn.rho.t[setting], gam.rho=gam.rho.t[setting],
                        nknots=nknots.t[setting])
    
    y[, , set, setting]       <- data$y
    r.inv.t.setting[, , set] <- data$r.inv
    knots.t.setting[, , , set] <- data$knots
    }
    
    r.inv.t[[setting]] <- r.inv.t.setting
    knots.t[[setting]] <- knots.t.setting
}

# remove simstudy "truth" settings to help diagnose errors.
save(y, r.inv.t, knots.t, ns, nt, s, nsets, 
	 x, ntest,		# covariate data that should be the same for all datasets
	 file='simdata.RData')
rm(list=ls())

# reload required functions into workspace
load(file='./simdata.RData')
source('../../R/mcmc.R')
source('../../R/densities.R')
source('../../R/auxfunctions.R')

save.image(file='simdata.RData')
