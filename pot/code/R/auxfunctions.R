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
rpotspatial <- function(nt, s, x, beta, sigma, delta, rho, nu, alpha, nknots=1){
    
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
  	sigma.t   <- sigma[t]
  	# print(sigma.t)
  	z.knots.t <- abs(rnorm(n=nknots, mean=0, sd=sqrt(sigma.t)))
    z.sites.t <- ZBySites(z.knots.t, partition[, t])
    
    x.beta.t <- x[, t, ] %*% beta
    
    # nugget <- (1 - alpha) * (sigma.t * (1 - delta^2))
    # p.sill <- alpha * (sigma.t * (1 - delta^2))
    
    nugget <- 1
    p.sill <- 1
    
    if (alpha != 0) {
      # cor    <- varcov.spatial(coords=s, cov.model="matern", nugget=nugget, 
      #                          cov.pars=c(p.sill, rho), kappa=nu)$varcov
      cor    <- alpha * matern(u=d, phi=rho, kappa=nu)
      v.spat <- t(chol(cor)) %*% rnorm(ns, 0, 1)
    } else {
      v.spat <- rep(0, ns)
    }
    
    # print(diag(cor))
    
    v.nug <- rnorm(ns, 0, sqrt(1 - alpha))
    v.t <- sqrt(sigma.t) * sqrt(1 - delta^2) * ( v.nug + v.spat) 
    # v.t    <- sqrt(sigma.t) * sqrt(1 - delta^2) * t(chol(cor)) %*% rnorm(ns)
    # v.t     <- t(chol(cor)) %*% rnorm(ns, mean=0, sd=sqrt(4))
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
#	cov(bool): do we need a variance term
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
#   z(nknots): vector of random effects with one entry per knot
#   partition(ns): partition membership vector
#
# Returns:
#   z.sites(ns): vector of random effects with one entry per site
################################################################
ZBySites <- function(z, partition){
  nknots  <- length(z)
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
#   prec.beta(p, p): prior precision of beta
#   e.beta(p): prior mean of beta
#   x(ns, nt, p): covariate array
#   y(ns, nt): observed data matrix
#   z(ns, nt): matrix of random effects
#   prec(ns, ns): precision matrix
#   delta(1): skewness parameter
#   sigma(nt): vector of daily variance
#   nt(1): number of days
#
# Returns:
#   list: 
#     vvv(p, p): posterior variance of beta
#     mmm(p): posterior mean of beta
################################################################
BetaPosterior <- function(prec.beta, e.beta, x, y, z,
						  prec, delta, sigma, nt){
  
  ns  <- nrow(y)
  p   <- length(e.beta)
  vvv <- prec.beta
  mmm <- e.beta
  
  for(t in 1:nt){
      x.t <- matrix(x[, t, ], ns, p)
      ttt  <- t(x.t) %*% prec / (sigma[t] * (1 - delta^2))
      vvv  <- vvv + ttt %*% x.t
      mmm  <- mmm + ttt %*% (y[, t] + delta * z[, t])
  }
	
  vvv <- chol2inv(chol(vvv))

  results <- list(vvv=vvv, mmm=mmm)

  return(results)
}

################################################################
# Arguments:
#   y(ns, nt): observed data matrix
#   x.beta(ns, nt): XBeta matrix
#   sigma(nt): vector of daily variance
#   delta(1): skewness parameter
#   prec(ns, ns): precision matrix
#   log.det(1): logdet(prec)
#   z.sites(ns, nt): matrix of random effects
#
# Returns:
#   llike(nt): (log)likelihood
################################################################
LLike <- function(y, x.beta, sigma, delta, prec, log.det, z.sites, log=TRUE){
  
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
    sigma.t <- sigma[t]
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



# ################################################################
# # Arguments:
# #   y.by.knots[nknots]: list of observed data matrices(ns.k, nt)
# #   mu.y[nknots]: list of E(Y) (ns.k, nt) by partition
# #	prec[nknots]: list of precision matrices(ns.k, ns.k)
# #   partition(ns): vector of partition membership
# #   nt(1): number of days
# #	nknots(1): number of knots
# #
# # Returns:
# #   ss(nknots, nt): sum of squares for each partition per day
# ################################################################
# SumSquares <- function(y.by.knots, mu.y, prec, partition, 
                       # nt, nknots){

  # ss <- matrix(NA, nknots, nt)

  # for (k in 1:nknots) {
    # these  <- which(partition == k) # identify sites to includes
    # ns.k   <- length(these)
    # if (ns.k > 0) {
      # res.k  <- y.by.knots[[k]] - mu.y[[k]]
      # prec.k <- prec[[k]]
      # for (t in 1:nt) {
        # ss[k, t] <- t(res.k[, t]) %*% prec.k %*% res.k[, t]
      # }
    # } else {
      # ss[k, ] <- rep(0, nt)
    # }
  # }
    
  # return(ss)
# }

# # ################################################################
# # Arguments:
# #   preds(yp, nt, iters): mcmc predictions at validation
# #                         locations
# #   probs(nprobs): sample quantiles for scoring
# #   validate(np, nt): validation data
# #
# # Returns:
# #   score(nprobs): a single quantile score per quantile
# ################################################################
# QuantScore <- function(preds, probs, validate){
  # nt <- ncol(validate)
  # np <- nrow(validate)
  # nprobs <- length(probs)
        
  # # apply gives nprobs x nsites. looking to find each site's quantile over all
  # # of the days.
  # pred.quants <- apply(preds, 1, quantile, probs=probs, na.rm=T)
    
  # scores.sites <- array(NA, dim=c(nprobs, np, nt))
    
  # for (q in 1:nprobs) {
    # diff <- pred.quants[q, ] - validate
    # i <- ifelse(diff >= 0, 1, 0)
    # scores.sites[q, , ] <- 2 * (i - probs[q]) * diff
  # }
    
  # scores <- apply(scores.sites, 1, mean, na.rm=T)

  # return(scores)
# }

# ################################################################
# # Arguments:
# #   preds(yp, nt, iters): mcmc predictions at validation
# #                         locations
# #   probs(nthreshs): sample quantiles for scoring
# #   validate(np, nt): validation data
# #
# # Returns:
# #   list:
# #     scores(nthreshs): a single brier score per threshold
# #     threshs(nthreshs): sample quantiles from dataset
# ################################################################
# BrierScore <- function(preds, probs, validate){
  # nthreshs <- length(probs)
  # thresholds <- quantile(validate, probs=probs, na.rm=T)
    
  # scores <- rep(NA, nthreshs)
  # for (b in 1:nthreshs) {
    # pat <- apply((preds > thresholds[b]), c(1, 2), mean)
    # ind <- validate < thresholds[b]
    # scores[b] <- mean((ind - pat)^2, na.rm=T)
  # }
    
  # return(scores)
# }

