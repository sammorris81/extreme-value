library(fields)
library(geoR)
library(mvtnorm)
library(evd)
#### begin setup file
rm(list=ls())
#### Necessary functions
source('./mcmc.R')
source('./densities.R')
source('./auxfunctions.R')

#########################################################################
# Simulation Study Settings
#########################################################################
#
# Settings:
#   alpha = 0.9
#   thresh = 0.9 Sample quantiles
#
#   Testing performance of gamma and fixed methods
#   nknots = c(1, 5) does partitioning help
#
#########################################################################

beta.y.t <- c(15, 2, -2)
sig2.y.t <- 1
rho.y.t <- 2
nu.y.t <- 0.5
alpha.y.t <- 0.9

xi.r.t <- 1     # first parameter for IG
sig.r.t <- 1    # second paratmeter for IG

nt <- 50
thresh <- 0

nknots <- c(1, 2)

#### Covariate data
s1 <- 1:7                       # Location x
s2 <- 1:7                       # Location y
S <- expand.grid(s1, s2)
ns <- nrow(S)  
X <- array(1,c(ns,nt,3))
for(t in 1:nt){
    X[,t,2]<-S[,1]
    X[,t,3]<-S[,2]
}

##############################################
# Dataset storage: 
#   settings: 1 knot or 3 knots
#   analysis: fixed vs gamma
#             1 knot vs 3 knots
#
#   
##############################################

nsets <- 1
nsettings <- 2 + 2 # fixed vs gamma, 1 knot vs 3 knots (for analysis settings)
nsettings <- nsettings + 1 #debug for setting all r = 1
threshold.1 <- threshold.3 <- full.1 <- full.3 <- array(NA, dim=c(ns, nt, nsets, nsettings))
thresh.gen.1 <- thresh.gen.3 <- matrix(NA, nrow=nsets, ncol=nsettings)
r.gen.1 <- array(NA, dim=c(nknots[1], nt, nsets, nsettings))
knots.gen.1 <- array(NA, dim=c(nknots[1], 2, nsets, nsettings))
r.gen.3 <- array(NA, dim=c(nknots[2], nt, nsets, nsettings))
knots.gen.3 <- array(NA, dim=c(nknots[2], 2, nsets, nsettings))
fixr <- c(F, F, F, F, T)
ntest <- floor(ns/2)
iters <- 20000
burn <- 10000

for(set in 1:nsets){
    for(setting in 1:nsettings){
    data <- rpotspatial(nt, S, X=X, beta.y=beta.y.t, rho.y=rho.y.t, nu.y=nu.y.t, 
                        alpha.y=alpha.y.t, thresh=thresh, xi.r=xi.r.t, 
                        sig.r=sig.r.t, fixr=fixr[setting], nknots=nknots[1])
    
    threshold.1[,,set,setting] <- data$y
    full.1[,,set,setting] <- data$y.full
    thresh.gen.1[set,setting] <- data$thresh.gen
    r.gen.1[,,set,setting] <- data$r
    knots.gen.1[,,set,setting] <- data$knots
                      
    data <- rpotspatial(nt, S, X=X, beta.y=beta.y.t, rho.y=rho.y.t, nu.y=nu.y.t, 
                        alpha.y=alpha.y.t, thresh=thresh, xi.r=xi.r.t, 
                        sig.r=sig.r.t, fixr=fixr[setting], nknots=nknots[2])
    
    threshold.3[,,set,setting] <- data$y
    full.3[,,set,setting] <- data$y.full
    thresh.gen.3[set,setting] <- data$thresh.gen
    r.gen.3[,,set,setting] <- data$r
    knots.gen.3[,,set,setting] <- data$knots
    
    }
}

# Remove simstudy "truth" settings to help diagnose errors.
save(threshold.1, full.1, thresh.gen.1, r.gen.1, knots.gen.1,
	 threshold.3, full.3, thresh.gen.3, r.gen.3, knots.gen.3, 
	 S, nsets, X,		# covariate data that should be the same for all datasets
	 file='simdata.RData')
rm(list=ls())

load(file='./simdata.RData')
source('./mcmc.R')
source('./densities.R')
source('./auxfunctions.R')

save.image(file='simdata.RData')

#### end setup file

#### begin pot file for 1 knots
rm(list=ls())
library(fields)
library(geoR)
library(mvtnorm)
library(evd)

load(file='./simdata.RData')
#source('./mcmc.R')
#source('./densities.R')
#source('./auxfunctions.R')

#### User specified data 

nknots <- 1						# how many knots to use for partition
thresh <- 0.0						# looking at predictions over this threshold
r.model <- "gamma"				# which model for random effect: gamma or fixed
iters <- 20000
burn <- 10000

#### Response data
setting <- 1 					# match with setting in manuscript
d <- 1							# will be assigned in for loop
y <- threshold.1[,,d,setting]
nt <- ncol(y)
ns <- nrow(y)
nsets <- dim(threshold.1)[3]
nsets <- 1 						# set to override the data

#### Cross-validation
y.full <- full.1[,,d,setting]
obs <- rep(c(T,F),ns)[1:ns]
y.o <- y[obs,]
y.full.o <- y.full[obs,]
X.o <- X[obs,,]
S.o <- S[obs,]
ntest <- nrow(y.full) - nrow(y.o)

Y.validate <- array(NA, dim=c(ntest, nt, nsets))	# validation data
Y.validate[,,d] <- y.p <- y[!obs,]
y.full.p <- y.full[!obs,]
X.p <- X[!obs,,]
S.p <- S[!obs,]

#### Sim study specific
pred <- array(NA, dim=c(ntest, nt, iters-burn, nsets))
params <- array(NA, dim=c(iters, 6, nsets))
beta <- array(NA, dim=c(iters-burn, 3, nsets))

outputfile <- paste("pot", setting, ".RData", sep="")
set.seed(setting)

start <- proc.time()	

set.seed(setting*100 + d)

# True settings for sim study
r.t <- rbind(r.gen.1[,,,setting])			# random effects
knots.t <- rbind(knots.gen.1[,,,setting])	# true location of knots
beta.y.t <- c(10, 0, 0)		# lm for E(Y)
sig2.y.t <- 1				# Fixed at 1 for identifiability
rho.y.t <- 2				# spatial range
nu.y.t <- 0.5				# matern smoothness
alpha.y.t <- 0.9			# spatial correlation strength (0 = ind, 1 = dep)
xi.r.t <- 1     			# first parameter for IG
sig.r.t <- 1    			# second paratmeter for IG


fit <- mcmc(y=y.o, S=S.o, X=X.o, Sp=S.p, Xp=X.p, thresh=thresh, r.model=r.model, 
			nknots=nknots, iters=iters, burn=burn, update=1000, thin=1, scale=T,
			debug=F,
			fixbeta=F, beta.init=beta.y.t,
			fixr=F, r.init = r.t,
			fixalpha=F, alpha.y.init = alpha.y.t,
			fixrhonu=F, logrho.y.init = log(rho.y.t), lognu.y.init = log(nu.y.t),
			fixxir=F, xi.r.init = xi.r.t,
			fixsigr=F, logsig.r.init = log(sig.r.t),
			fixknots=F, knots.init = knots.t
			)
#### end 1 knot


y=y.o; S=S.o; X=X.o; Sp=S.p; Xp=X.p; thresh=thresh; r.model=r.model; 
nknots=nknots; iters=iters; burn=burn; update=1000; thin=1; 

beta.y.m=0; beta.y.s=10; beta.init=0;
logrho.y.m=log(0.5); logrho.y.s=1; logrho.y.init=log(0.5)
lognu.y.m=log(0.5); lognu.y.s=0.01; lognu.y.init=log(0.5)
alpha.y.a=1; alpha.y.b=1; alpha.y.init=0.5    
mu.r.m=0; mu.r.s=10; mu.r.init=0
logsig.r.m=0; logsig.r.s=1; logsig.r.init=0
xi.r.m=0; xi.r.s=0.2; xi.r.init=0
iters=iters; burn=burn; update=100; thin=1
debug=T
scale = T;
fixbeta = F; fixalpha = F; fixrhonu = F; 
fixr = T; r.init = r.t; fixxir = T; fixsigr = T; fixknots= F; 

#### begin 3 knots

#### begin pot file for 3 knots
rm(list=ls())
library(fields)
library(geoR)
library(mvtnorm)
library(evd)

load(file='./simdata.RData')
#source('./mcmc.R')
#source('./densities.R')
#source('./auxfunctions.R')

#### User specified data 

nknots <- 12						# how many knots to use for partition
thresh <- 0.0					# looking at predictions over this threshold
r.model <- "gamma"				# which model for random effect: gamma or fixed
iters <- 20000
burn <- 10000

#### Response data
setting <- 1 					# match with setting in manuscript
d <- 1							# will be assigned in for loop
y <- threshold.1[,,d,setting]
nt <- ncol(y)
ns <- nrow(y)
nsets <- dim(threshold.3)[3]
nsets <- 1						# set to override the data

#### Cross-validation
y.full <- full.1[,,d,setting]
obs <- rep(c(T,F),ns)[1:ns]
y.o <- y[obs,]
y.full.o <- y.full[obs,]
X.o <- X[obs,,]
S.o <- S[obs,]
ntest <- nrow(y.full) - nrow(y.o)

Y.validate <- array(NA, dim=c(ntest, nt, nsets))	# validation data
Y.validate[,,d] <- y.p <- y[!obs,]
y.full.p <- y.full[!obs,]
X.p <- X[!obs,,]
S.p <- S[!obs,]

#### Sim study specific
pred <- array(NA, dim=c(ntest, nt, iters-burn, nsets))
params <- array(NA, dim=c(iters, 7, nsets))
beta <- array(NA, dim=c(iters-burn, 3, nsets))

outputfile <- paste("pot", setting, ".RData", sep="")
set.seed(setting)

start <- proc.time()	

set.seed(setting*100 + d)

# True settings for sim study
r.t <- rbind(r.gen.3[,,,setting])			# random effects
knots.t <- rbind(knots.gen.3[,,,setting])	# true location of knots
beta.y.t <- c(10, 0, 0)		# lm for E(Y)
sig2.y.t <- 1				# Fixed at 1 for identifiability
rho.y.t <- 2				# spatial range
nu.y.t <- 0.5				# matern smoothness
alpha.y.t <- 0.9			# spatial correlation strength (0 = ind, 1 = dep)
xi.r.t <- 1     			# first parameter for IG
sig.r.t <- 1    			# second paratmeter for IG


fit <- mcmc(y=y.o, S=S.o, X=X.o, Sp=S.p, Xp=X.p, thresh=thresh, r.model=r.model, 
			nknots=nknots, iters=iters, burn=burn, update=1000, thin=1, 
			debug=F,
			fixbeta=F, beta.init=beta.y.t,
			fixr=F, r.init = r.t,
			fixalpha=F, alpha.y.init = alpha.y.t,
			fixrhonu=F, logrho.y.init = log(rho.y.t), lognu.y.init = log(nu.y.t),
			fixxir=F, xi.r.init = xi.r.t,
			fixsigr=F, logsig.r.init = log(sig.r.t),
			fixknots=F, knots.init = knots.t
			)
			
#### end 3 knots

y=y.o; S=S.o; X=X.o; Sp=S.p; Xp=X.p; thresh=thresh; r.model=r.model; 
nknots=nknots; iters=iters; burn=burn; update=1000; thin=1; 

beta.y.m=0; beta.y.s=10; beta.init=0;
logrho.y.m=log(1); logrho.y.s=1; logrho.y.init=log(1)
lognu.y.m=log(0.5); lognu.y.s=0.01; lognu.y.init=log(0.5)
alpha.y.a=1; alpha.y.b=1; alpha.y.init=0.5    
mu.r.m=0; mu.r.s=10; mu.r.init=0
logsig.r.m=0; logsig.r.s=1; logsig.r.init=0
xi.r.m=0; xi.r.s=0.2; xi.r.init=0
iters=iters; burn=burn; update=100; thin=1
debug=T; scale=T;
fixbeta = F; fixxir = F; fixalpha = F; fixrhonu = F; fixr = F;
fixsigr = F; fixknots= F; 


#### Make sure the partition stuff gives back same as non-partitioned

beta.y.t <- c(10, 0, 0)
sig2.y.t <- 1
rho.y.t <- 2
nu.y.t <- 0.5
alpha.y.t <- 0.9

xi.r.t <- 1     # first parameter for IG
sig.r.t <- 1    # second paratmeter for IG

nt <- 10
thresh <- 0.9

nknots <- c(1, 2)

#### Covariate data
s1 <- 1:3                       # Location x
s2 <- 1:4                       # Location y
S <- expand.grid(s1, s2)
S.scale <- S
ns <- nrow(S)  
X <- array(1,c(ns,nt,3))
for(t in 1:nt){
    X[,t,2]<-S[,1]
    X[,t,3]<-S[,2]
}


# Initial setup
ns <- nrow(S)
d <- rdist(S)
diag(d) <- 0
min1 <- min(S[,1])
max1 <- max(S[,1])
min2 <- min(S[,2])
max2 <- max(S[,2])

# generate the knot locations
knots1 <- matrix(NA, nknots[1], 2)
knots1[,1] <- runif(nknots[1], min1, max1)
knots1[,2] <- runif(nknots[1], min2, max2)

g1 <- membership(S, knots1)


# generate the knot locations
knots2 <- matrix(c(1.25, 1.25, 2.75,2.75), nknots[2], 2, byrow=T)

g2 <- membership(S, knots2)

# sort the locations by partition
which(g2$partition==2)
tempS <- S[6,]
S[6,] <- S[7,]
S[7,] <- tempS

d <- rdist(S)
diag(d) <- 0

# Generate the spatial correlation matrix.
Sig1 <- spatcor(d, alpha.y.t, rho.y.t, nu.y.t, g1$blockcor, g1$partition)
SIG1 <- Sig1$SIG
PREC1 <- Sig1$PREC


g2 <- membership(S, knots2)

# Generate the spatial correlation matrix.
Sig2 <- spatcor(d, alpha.y.t, rho.y.t, nu.y.t, g2$blockcor, g2$partition)
SIG2 <- Sig2$SIG
PREC2 <- Sig2$PREC

round(SIG1, 3)
round(SIG2, 3)
round(PREC1, 3)
round(PREC2, 3)










