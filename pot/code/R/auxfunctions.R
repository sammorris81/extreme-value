####################################################################
# Auxiliary Functions
####################################################################

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
  
  lower.u <- ifelse( mn / sd > 5 & lower == 0, 0, lower.u )
  U <- runif(length(mn), lower.u, upper.u)
  y <- qnorm(U, mn, sd)
  
  return(y)
}

#########################################################################
# Arguments:
#   nt(1): number of days
#   s(ns, 2): spatial locations
#   x(ns, nt, np): Matrix of spatial covariates
#   beta.y(np): parameters for mu.y
#   sigma.sites(ns, nt): matrix of daily variance at each site
#   delta(1): skewness parameter
#   rho(1): spatial range parameter
#   nu.y(1): smoothness parameter for Matern covariance
#   alpha.y(1): controls proportion of spatial to non-spatial covariance
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
rpotspatial <- function(nt, s, x, beta, sigma.sites, delta, rho, nu, alpha, nknots=1){
    
  # initial setup
  ns      <- nrow(s)
  y       <- matrix(NA, ns, nt)
  z.knots <- matrix(NA, nrow=nknots, ncol=nt)
  z.sites <- matrix(NA, nrow=ns, ncol=nt)
  v       <- matrix(NA, nrow=ns, ncol=nt)
  knots   <- vector(mode="list", length=nt)
  
  d       <- rdist(s)
  diag(d) <- 0
  min1    <- min(s[, 1])
  max1    <- max(s[, 1])
  min2    <- min(s[, 2])
  max2    <- max(s[, 2])
  
  for (t in 1:nt) {
    # generate the knot locations
    knots[[t]]      <- matrix(NA, nknots, 2)
    knots[[t]][, 1] <- runif(nknots, min1, max1)
    knots[[t]][, 2] <- runif(nknots, min2, max2)
  }
  
  partition    <- Membership(s=s, knots=knots)
   
  for (t in 1:nt) {
  	sigma.t   <- sigma.sites[, t]
  	z.knots.t <- sqrt(sigma.t) * abs(rnorm(n=nknots, mean=0, sd=1))
    z.sites.t <- ZBySites(z.knots.t, partition[, t])
    
    x.beta.t <- x[, t, ] %*% beta
     
    if (alpha != 0) {
      cor    <- alpha * matern(u=d, phi=rho, kappa=nu)
      v.spat <- t(chol(cor)) %*% rnorm(ns, 0, 1)
    } else {
      v.spat <- rep(0, ns)
    }
    
    v.nug <- rnorm(ns, 0, sqrt(1 - alpha))
    v.t <- sqrt(sigma.t) * sqrt(1 - delta^2) * ( v.nug + v.spat) 
    y[, t]       <- x.beta.t + delta * z.sites.t + v.t
    z.knots[, t] <- z.knots.t
    z.sites[, t] <- z.sites.t
    v[, t]       <- v.t
    
  }
  
  results <- list(y=y, z.knots=z.knots, z.sites=z.sites, knots=knots, v=v, cor=cor)
    
  return(results)
}

################################################################
# Arguments:
#   s(ns, 2): spatial locations
#   knots[nt](nknots, 2): list of knot locations matrices
#
# Returns: 
#   partition(ns, nt): partition membership matrix
################################################################
Membership <- function(s, knots){
  ns <- nrow(s)
  nt <- length(knots)
  nknots <- nrow(knots[[1]])
  partition <- matrix(NA, ns, nt) # a number letting us know which partition.
  
  # membership matrix
  for (t in 1:nt) {
    d <- rdist(s, knots[[t]])
    partition[, t] <- apply(d, 1, which.min)
  }
  
  return(partition)
}

################################################################
# Arguments:
#   sigma.knots(nknots, nt): matrix of daily variance at each knot
#   partition(ns, nt): partition membership matrix
#   nknots(1): the number of knots/partitions to have over space
#
# Returns: 
#   sigma.sites(ns, nt): matrix of daily variance at each site
################################################################
SigmaSites <- function(sigma.knots, partition, nknots){
  ns <- nrow(partition)
  nt <- ncol(partition)
  
  sigma.sites <- matrix(NA, ns, nt)
  
  for (t in 1:nt) {
    partition.t <- partition[, t]
    for (k in 1:nknots) {
      these <- which(partition.t == k)
      sigma.sites[these, t] <- sigma.knots[k, t]
    }
  }
  
  return(sigma.sites)
}


################################################################
# Arguments:
#   s(ns, 2): spatial locations
#
# Returns:
#   s.scale(ns, 2): locations scaled to be in [0, 1] x [0, 1]
################################################################
ScaleLocs <- function(s){
  x.min <- min(s[, 1]); x.max <- max(s[, 1]); x.range <- x.max - x.min
  y.min <- min(s[, 2]); y.max <- max(s[, 2]); y.range <- y.max - y.min

  s.x <- (s[,1] - x.min) / x.range
  s.y <- (s[,2] - y.min) / y.range
	
  s.scale <- cbind(s.x, s.y)
	
  return(s.scale)   
}



################################################################
# Arguments:
#   d(ns, ns): distance between observations
#   alpha(1): controls proportion of spatial to non-spatial 
#             covariance (0: ind, 1: high spatial corr)
#   rho(1): spatial range
#   nu(1): matern smoothness parameter
#   eps(1): small amount for numerical stability
#
# Returns:
#   list: 
#     prec(ns, ns): precision matrix
#     log.det(nknots): logdet(prec)
#     sig(ns, ns): correlation matrix
################################################################
SpatCor <- function(d, alpha, rho, nu=0.5, eps=10^(-5)){
  q <- sig <- list()

  sig       <- CorFx(d, alpha, rho, nu, cov=FALSE)
  sig.chol  <- chol(sig)
  diag.chol <- ifelse(diag(sig.chol) < eps, eps, diag(sig.chol))
  log.det   <- -2 * sum(log(diag.chol))
  prec      <- chol2inv(sig.chol)
  
  results <- list(prec=prec, log.det=log.det, sig=sig)

  return(results)

}

################################################################
# Arguments:
#   d(ns, ns): distance between observations
#   alpha(1): controls proportion of spatial to non-spatial 
#             covariance (0: ind, 1: high spatial corr)
#   rho(1): spatial range
#   nu(1): matern smoothness parameter
#	  cov(bool): do we need a variance term
#
# Returns:
#   cor(ns, nt): matern correlation
################################################################
CorFx <- function(d, alpha, rho, nu, cov=F){
    
  ns <- nrow(d)
  cor <- alpha * matern(d, rho, nu)
  if (!cov) { # S12 doesn't need diagonal variance term
    cor <- cor + (1 - alpha) * diag(rep(1, ns)) 
  } 
    
  return(cor)
}

################################################################
# Arguments:
#   z.knots(nknots): vector of random effects with one entry per knot
#   partition(ns): partition membership vector
#
# Returns:
#   z.sites(ns): vector of random effects with one entry per site
################################################################
ZBySites <- function(z.knots, partition){
  nknots  <- length(z.knots)
  ns      <- length(partition)
  z.sites <- rep(NA, ns)
  
  for (k in 1:nknots) {
  	these <- which(partition == k)
  	z.sites[these] <- z[k]
  }
  
  return(z.sites)
}

################################################################
# Arguments:
#   res(ns, nt)
#   prec(ns, ns): precision matrix
#   sigma.sites(ns, nt): matrix of daily variance at each site
#
# Returns:
#   ss(nt): vector of daily sums of squares
################################################################
SumSquares <- function(res, prec, sigma.sites){
  nt <- ncol(res)
  ss <- rep(NA, nt)
  
  for (t in 1:nt) {
    sigma.sites.t <- sqrt(sigma.sites[, t])
    ss[t] <- t(res[, t] / sigma.sites.t) %*% prec %*% res[, t] / sigma.sites.t
  }
  
  return(ss)
}

################################################################
# Arguments:
#   y(ns, nt): observed data matrix
#   x.beta(ns, nt): XBeta matrix
#   sigma.sites(ns, nt): matrix of daily variance at each site
#   delta(1): skewness parameter
#   prec(ns, ns): precision matrix
#   log.det(1): logdet(prec)
#   z.sites(ns, nt): matrix of random effects
#
# Returns:
#   llike(nt): (log)likelihood
################################################################
LLike <- function(y, x.beta, sigma.sites, delta, prec, log.det, z.sites, log=TRUE){
  
  if (missing(y)) {
    stop("y must be defined")
  } 
  ns <- nrow(y)
  nt <- ncol(y)
  if (missing(x.beta)) {
    stop("x.beta must be defined")
  } else if (nrow(x.beta) != ns || ncol(x.beta) != nt) {
    stop("x.beta and y must have conforming size")
  }
  
  if (missing(sigma)) {
    stop("sigma must be defined")
  } else if (length(sigma) != nt) {
    stop("need a sigma for each day")
  }
  
  if (missing(delta)) {
    stop("delta must be defined")
  }
  
  if (missing(log.det)) {
    stop("log.det must be defined")
  }
      
  log.like  <- rep(NA, nt)

  for (t in 1:nt) {
    sigma.t <- sigma.sites[, t]
    y.t <- y[, t]
    mu.t <- x.beta[, t] + delta * z.sites[, t]
    ss.t <- t(y.t - mu.t) %*% prec %*% (y.t - mu.t)
    log.like[t] <- -0.5 * ns * (log(sigma.t) + log(1 - delta^2)) + 0.5 * log.det -
                    0.5 * ss.t / (sigma.t * (1 - delta^2))
  }
    
  if(!log){
    log.like <- exp(log.like)
  }
    
  return(log.like)
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
