###############################################################################
#### Basic setup                                                           ####
###############################################################################
#setwd('~/svn-repos/dissertation/extreme-value/code')
#setwd('S:/Documents/Dissertation/extreme-value/code')

library(fields)
library(emulator)
library(evd)
library(POT)
library(geoR)

#### mcmc.R - MCMC algorithm
source('./mcmc.R')

#### densities.R - Required density calculations for mcmc
source('./densities.R')

#### Spatial GP code
source('./spatial_gp.R')

###############################################################################
#### Data settings                                                         #### 
###############################################################################

#### Model parameters                

## Latent variable
alpha <- c(.3, .5, .7)  # Spatial correlation 
                      # (0 = spatial dep, 1 = spatial ind)     
logrange <- 0       # Kernel bandwidth 

## Regression
xi <- .1
nt <- 10       				    # Number of days 

#### Two options for beta terms
beta.1 <- c(3.458, 0, 0)  # Intercept, Lat, Long 5%
beta.2 <- c(1.993, 0, 0)  # Intercept, Lat, Long 15%

#### Three different grid sizes for checking timing
s1 <- 1:6       				# Location x
s2 <- 1:6      				    # Location y
S.1 <- expand.grid(s1, s2)
knots.1 <- S.1      			# Knots at grid points  
n.1 <- nrow(S.1)  
X.hm.1 <- as.matrix(cbind(1, S.1), n.1, 3)

s1 <- 1:10       				# Location x
s2 <- 1:10      			    # Location y
S.2 <- expand.grid(s1, s2)
knots.2 <- S.2      			# Knots at grid points               
n.2 <- nrow(S.2)
X.hm.2 <- as.matrix(cbind(1, S.2), n.2, 3)

s1 <- 1:15       				# Location x
s2 <- 1:15     				    # Location y
S.3 <- expand.grid(s1, s2)
knots.3 <- S.3      			# Knots at grid points               
n.3 <- nrow(S.3)
X.hm.3 <- as.matrix(cbind(1, S.3), n.3, 3)

#### Thresholds
thresholds.1.1 <- (1 + xi * X.hm.1 %*% beta.1)^(1/xi)
thresholds.1.2 <- (1 + xi * X.hm.1 %*% beta.2)^(1/xi)
thresholds.2.1 <- (1 + xi * X.hm.2 %*% beta.1)^(1/xi)
thresholds.2.2 <- (1 + xi * X.hm.2 %*% beta.2)^(1/xi)
thresholds.3.1 <- (1 + xi * X.hm.3 %*% beta.1)^(1/xi)


#### MCMC Settings
iters <- 25000
burn  <- 15000
#ntest <- 0
#ntrain <- nrow(S) - ntest

nsets <- 1

y.1.1 <- y.1.2 <- latent.1.1 <- latent.1.2 <- array(data=NA, dim=c(nrow(S.1), nt, (nsets)))
y.2.1 <- y.2.2 <- latent.2.1 <- latent.2.2 <- array(data=NA, dim=c(nrow(S.2), nt, (nsets)))
y.3.1 <- latent.3.1 <- array(data=NA, dim=c(nrow(S.3), nt, (nsets)))

set.seed(2087)

for(i in 1:(nsets)){
	y.gen <- rgevcond(
		nreps.f=nt, S.f=S.1, knots.f=knots.1, alpha.f=alpha[2], logrange.f=logrange, theta.f=NULL
	) 
	y.1.1[,,i] <- as.numeric(y.gen$y > rep(thresholds.1.1, nt))
	latent.1.1[,,i] <- y.gen$y

	y.gen <- rgevcond(
		nreps.f=nt, S.f=S.1, knots.f=knots.1, alpha.f=alpha[2], logrange.f=logrange, theta.f=NULL
	) 
	y.1.2[,,i] <- as.numeric(y.gen$y > rep(thresholds.1.2, nt))
	latent.1.2[,,i] <- y.gen$y

	y.gen <- rgevcond(
		nreps.f=nt, S.f=S.2, knots.f=knots.2, alpha.f=alpha[2], logrange.f=logrange, theta.f=NULL
	) 
	y.2.1[,,i] <- as.numeric(y.gen$y > rep(thresholds.2.1, nt))
	latent.2.1[,,i] <- y.gen$y
	
	y.gen <- rgevcond(
		nreps.f=nt, S.f=S.2, knots.f=knots.2, alpha.f=alpha[2], logrange.f=logrange, theta.f=NULL
	) 
	y.2.2[,,i] <- as.numeric(y.gen$y > rep(thresholds.2.2, nt))
	latent.2.2[,,i] <- y.gen$y
	
	y.gen <- rgevcond(
		nreps.f=nt, S.f=S.3, knots.f=knots.3, alpha.f=alpha[2], logrange.f=logrange, theta.f=NULL
	) 
	y.3.1[,,i] <- as.numeric(y.gen$y > rep(thresholds.3.1, nt))
	latent.3.1[,,i] <- y.gen$y
}

y.1 <- cbind(y.1.1[,7,], y.1.2[,2,])
latent.1 <- cbind(latent.1.1[,7,], latent.1.2[,2,])
y.2 <- cbind(y.2.1[,5,], y.2.2[,1,])
latent.2 <- cbind(latent.2.1[,5,], latent.2.2[,1,])
y.3 <- cbind(y.3.1[,3,])
latent.3 <- cbind(latent.3.1[,3,])

#rm(alpha, logrange, xi, beta)
#save.image(file="simdata.RData")
save.list <- as.list(y.1, latent.1, y.2, latent.2, y.3, latent.3)
save(list=save.list, file="timetrial.RData")