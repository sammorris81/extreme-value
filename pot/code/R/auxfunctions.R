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
#   tau.alpha.t(1): hyperparameter for tau terms
#   tau.beta.t(1): hyperparameter for tau terms
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
rpotspatial <- function(nt, s, x, beta, tau.alpha.t, tau.beta.t, 
                        delta, rho, nu, alpha, nknots=1){
    
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
  tau.knots <- rgamma(nknots * nt, tau.alpha.t, tau.beta.t)
  sigma.knots <- 1 / tau.knots
  z.knots <- sqrt(sigma.knots) * abs(rnorm(n=(nknots * nt), mean=0, sd=1))
  
  tau.knots <- matrix(tau.knots, nrow=nknots, ncol=nt)
  tau.sites <- TauSites(tau.knots, partition, nknots)
  sigma.sites <- 1 / tau.sites
  z.knots <- matrix(z.knots, nrow=nknots, ncol=nt)
  z.sites <- ZBySites(z.knots, partition, nknots)
  
  for (t in 1:nt) {
    x.beta.t <- x[, t, ] %*% beta
    
    sigma.sites.t <- sigma.sites[, t]
    cor.mtx       <- SpatCor(d=d, alpha=alpha, rho=rho, nu=nu)
    cor           <- cor.mtx$sig
    cov.t         <- diag(sqrt(sigma.sites.t)) %*% cor %*% diag(sqrt(sigma.sites.t))
     
    # if (alpha != 0) {
      # cor    <- alpha * matern(u=d, phi=rho, kappa=nu)
      # v.spat <- t(chol(cor)) %*% rnorm(ns, 0, 1)
    # } else {
      # v.spat <- rep(0, ns)
    # }
    
    z.sites.t <- z.sites[, t]
    # v.nug     <- rnorm(ns, 0, sqrt(1 - alpha))
    v.t       <- sqrt(1 - delta^2) * t(chol(cov.t)) %*% rnorm(ns, 0, 1) 
    y[, t]    <- x.beta.t + delta * z.sites.t + v.t
    v[, t]    <- v.t
    
  }
  
  results <- list(y=y, z.knots=z.knots, z.sites=z.sites, knots=knots, v=v, cor=cor,
                  tau.knots=tau.knots, tau.sites=tau.sites)
    
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
#   tau.knots(nknots, nt): matrix of daily precision at each knot
#   partition(ns, nt): partition membership matrix
#   nknots(1): the number of knots/partitions to have over space
#
# Returns: 
#   tau.sites(ns, nt): matrix of daily precision at each site
################################################################
TauSites <- function(tau.knots, partition, nknots){
  ns <- nrow(partition)
  nt <- ncol(partition)
  
  tau.sites <- matrix(NA, ns, nt)
  
  for (t in 1:nt) {
    partition.t <- partition[, t]
    for (k in 1:nknots) {
      these <- which(partition.t == k)
      tau.sites[these, t] <- tau.knots[k, t]
    }
  }
  
  return(tau.sites)
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
#     prec.cor(ns, ns): precision matrix (correlation)
#     logdet.prec(nknots): logdet(prec)
#     sig(ns, ns): correlation matrix
################################################################
SpatCor <- function(d, alpha, rho, nu=0.5, eps=10^(-5)){
  q <- sig <- list()

  sig         <- CorFx(d, alpha, rho, nu, cov=FALSE)
  sig.chol    <- chol(sig)
  diag.chol   <- ifelse(diag(sig.chol) < eps, eps, diag(sig.chol))
  logdet.prec <- -2 * sum(log(diag.chol))
  prec.cor    <- chol2inv(sig.chol)
  
  results <- list(prec.cor=prec.cor, logdet.prec=logdet.prec, sig=sig)

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
#   z.knots(nknots, nt): matrix of random effects with one entry per knot/day
#   partition(ns, nt): partition membership matrix
#
# Returns:
#   z.sites(ns, nt): matrix of random effects with one entry per site/day
################################################################
ZBySites <- function(z.knots, partition, nknots){
  if (nrow(z.knots) != nknots) {
    stop("z.knots has the wrong number of rows")
  }
  
  if (ncol(z.knots) != ncol(partition)) {
    stop("z.knots and partition should have the same number of columns")
  }
  
  ns      <- nrow(partition)
  nt      <- ncol(partition)
  z.sites <- matrix(NA, ns, nt)
  
  for (t in 1:nt) {
    for (k in 1:nknots) {
  	  these <- which(partition[, t] == k)
  	  z.sites[these, t] <- z.knots[k, t]
    }
  }
  
  
  return(z.sites)
}

################################################################
# Arguments:
#   res(ns, nt): matrix of residuals
#   prec.cor(ns, ns): precision matrix (correlation)
#   tau.sites(ns, nt): matrix of daily precision at each site
#
# Returns:
#   ss(nt): vector of daily sums of squares
################################################################
SumSquares <- function(res, prec.cor, tau.sites){
  if ((nrow(res) != nrow(tau.sites)) | (ncol(res) != ncol(tau.sites))) {
    stop("res and tau.sites do not have the same dimensions")
  }
  
  nt <- ncol(res)
  ss <- rep(NA, nt)
  
  for (t in 1:nt) {
    tau.sites.t <- sqrt(tau.sites[, t])
    ss[t] <- t(res[, t] * tau.sites.t) %*% prec.cor %*% (res[, t] * tau.sites.t)
  }
  
  return(ss)
}

################################################################
# Arguments:
#   y(ns, nt): observed data matrix
#   x.beta(ns, nt): XBeta matrix
#   tau.sites(ns, nt): matrix of daily precision at each site
#   delta(1): skewness parameter
#   prec.cor(ns, ns): precision matrix (correlation)
#   logdet.prec(1): logdet(prec)
#   z.sites(ns, nt): matrix of random effects
#
# Returns:
#   llike(nt): (log)likelihood
################################################################
LLike <- function(y, x.beta, tau.sites, delta, prec.cor, logdet.prec, z.sites, log=T){
  
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
  
  if (missing(tau.sites)) {
    stop("tau.sites must be defined")
  } else if ((nrow(tau.sites) != ns) | (ncol(tau.sites) != nt)) {
    stop("need a tau for each day/knot")
  }
  
  if (missing(delta)) {
    stop("delta must be defined")
  }
  
  if (missing(logdet.prec)) {
    stop("logdet.prec must be defined")
  }
      
  log.like  <- rep(NA, nt)

  for (t in 1:nt) {
    tau.t <- sqrt(tau.sites[, t])
    y.t <- y[, t]
    mu.t <- x.beta[, t] + delta * z.sites[, t]
    std.res.t <- (y.t - mu.t) * tau.t
    ss.t <- t(std.res.t) %*% prec.cor %*% (std.res.t)
    log.like[t] <- 0.5 * (sum(log(tau.t)) - ns * (log(1 - delta^2))) + 
                   0.5 * logdet.prec - 0.5 * ss.t / (1 - delta^2)
  }
  
  if(!log){
    log.like <- exp(log.like)
  }
    
  return(log.like)
}

################################################################
# Arguments:
#   preds(yp, iters): mcmc predictions at validation
#                         locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
QuantScore <- function(preds, probs, validate){
  np <- length(validate)
  nprobs <- length(probs)
        
  # apply gives nprobs x n. looking to find each site's quantile over all
  # of the days.
  pred.quants <- apply(preds, 1, quantile, probs=probs, na.rm=T)
    
  scores.sites <- matrix(NA, nrow=nprobs, ncol=np)
    
  for (q in 1:nprobs) {
    diff <- pred.quants[q] - validate
    i <- ifelse(diff >= 0, 1, 0)
    scores.sites[q, ] <- 2 * (i - probs[q]) * diff
  }
    
  scores <- apply(scores.sites, 1, mean, na.rm=T)

  return(scores)
}

################################################################
# Arguments:
#   preds(yp, iters): mcmc predictions at validation
#                         locations
#   probs(nthreshs): sample quantiles for scoring
#   validate(np): validation data
#
# Returns:
#   list:
#     scores(nthreshs): a single brier score per threshold
################################################################
BrierScore <- function(preds, probs, validate){
  nthreshs <- length(probs)
  thresholds <- quantile(validate, probs=probs, na.rm=T)
    
  scores <- rep(NA, nthreshs)
  for (b in 1:nthreshs) {
    pat <- apply((preds > thresholds[b]), 1, mean)
    ind <- validate < thresholds[b]
    scores[b] <- mean((ind - pat)^2, na.rm=T)
  }
    
  return(scores)
}

