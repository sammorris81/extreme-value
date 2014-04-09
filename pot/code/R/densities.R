#########################################################################
# Arguments:
#   nt(1): number of days
#   s(ns, 2): spatial locations
#   x(ns, nt, np): Matrix of spatial covariates
#   beta.y(np): parameters for mu.y
#   nu.y(1): smoothness parameter for Matern covariance
#   alpha.y(1): controls proportion of spatial to non-spatial covariance
#   thresh(1): percentage of data above threshold in GPD
#   xi.r(1): a parameter for IG
#   sig.r(1): b parameter for IG
#	mixprob(1): 0: all IG random effects; 1: all MVN random effects
#	mvn.rho(1): range parameter for MVN correlation
#	gam.rho(1): range parameter for IG correlation
#   nknots(1): the number of knots/partitions to have over space
#
# Returns:
#   list:
#		y(ns, nt): data
#		r.inv(nknots, nt): random effects for each knot and day
#		thresh.gen(1): sample quantile at which data is thresholded
#		knots(nknots, 2): knot locations
#       fixr(nknots, nt): whether it was a gaussian or gamma random effect
#########################################################################
rpotspatial <- function(nt, s, x, beta.y, nu.y, alpha.y, thresh, 
						xi.r, sig.r, mixprob=0.5, 
                        mvn.rho, gam.rho, nknots=1){
    
  # initial setup
  ns      <- nrow(s)
  d       <- rdist(s)
  diag(d) <- 0
  min1    <- min(s[, 1])
  max1    <- max(s[, 1])
  min2    <- min(s[, 2])
  max2    <- max(s[, 2])
    
  # generate the knot locations
  knots      <- matrix(NA, nknots, 2)
  knots[, 1] <- runif(nknots, min1, max1)
  knots[, 2] <- runif(nknots, min2, max2)
  
  y.temp     <- matrix(0, ns, nt)  # only used so Membership behaves properly
  member     <- Membership(s=s, knots=knots, y=y.temp, x=x)
  partition  <- member$partition
  x.by.knots <- member$x.by.knots
  mu.y       <- ExpectY(partition=partition, x.by.knots=x.by.knots, 
                        beta.y=beta.y, nknots=nknots)
    
  y     <- matrix(NA, ns, nt)  # data storage
  r.inv <- matrix(NA, nrow=nknots, ncol=nt)
  fixr  <- matrix(NA, nrow=nknots, ncol=nt)
  
  for (k in 1:nknots) {
    these  <- which(partition == k)
    ns.k   <- length(these)
    if(ns.k == 1){
    	s.k <- matrix(s[these, ], 1, 2)
    } else {
    	s.k <- s[these, ]
    }
    mu.y.k <- matrix(mu.y[[k]], nrow=ns.k, ncol=nt)
    
    if (ns.k > 0) {
      for(t in 1:nt){  # random effects differ by day for each knot
        fixr.kt <- rbinom(1, 1, mixprob)
        if (fixr.kt) {  # gaussian
          r.inv.kt <- 2
          rho.kt   <- mvn.rho
        } else {
      	  r.inv.kt <- rgamma(1, xi.r, sig.r)
          rho.kt   <- gam.rho
        }
        nug.kt    <- (1 - alpha.y) / r.inv.kt
        p.sill.kt <- alpha.y / r.inv.kt
        mu.y.kt   <- mu.y.k[, t]
        sig.kt    <- varcov.spatial(s.k, cov.model="matern", 
                                    nugget=nug.kt, cov.pars=c(p.sill.kt, rho.kt),
                                    kappa=nu.y)$varcov
        y[these,t] <- rmvnorm(1, mean=mu.y.kt, sigma=sig.kt)
        r.inv[k, t] <- r.inv.kt
        fixr[k, t] <- fixr.kt
      }
    }
    
  }

  thresh.gen <- quantile(y, probs=thresh)
        
  results <- list(y=y, r.inv=r.inv, thresh.gen=thresh.gen, knots=knots, fixr=fixr)
  
  return(results)
}

################################################################
# IG(a, b) density function
################################################################
dInvG <- function(x, a, b, log=T){
  lll <- -(a + 1) * log(x) - b / x
  if (!log) {
    lll<-exp(lll)
  }
    
  return(lll)
}

#########################################################################
# Arguments:
#   mn(nt): mean
#   sd(nt, nt): standard deviation
#   lower(1): lower truncation point (default=-Inf) 
#   upper(1): upper truncation point (default=Inf)
#	fudge(1): small number for numerical stability (used for lower bound)
#
# Returns:
#   y(nt): truncated normal data
#########################################################################

rTNorm <- function(mn, sd, lower=-Inf, upper=Inf, fudge=0){
  lower.u <- pnorm(lower, mn, sd)
  upper.u <- pnorm(upper, mn, sd)

  U <- runif(length(mn), lower.u, upper.u)
  y <- qnorm(U, mn, sd)
    
  return(y)
}