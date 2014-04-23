####################################################################
# Auxiliary Functions
####################################################################

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
  log.det  <- rep(0, nknots)
  
  sig       <- CorFx(d, alpha, rho, nu, cov=FALSE)
  sig.chol  <- chol(sig)
  diag.chol <- ifelse(diag(sig.chol) < eps, eps, diag(sig.chol))
  log.det   <- 2 * sum(log(diag.chol))
  log.det   <- 1 / log.det
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
#   z(nknots): vector of random effects
#   partition(ns): partition membership vector
#
# Returns:
#   z.sites(ns): matrix of random effects for mu
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
#   y.by.knots[nknots]: list of observed data matrices(ns.k, nt)
#   mu.y[nknots]: list of E(Y) (ns.k, nt) by partition
#	prec[nknots]: list of precision matrices(ns.k, ns.k)
#   partition(ns): vector of partition membership
#   nt(1): number of days
#	nknots(1): number of knots
#
# Returns:
#   ss(nknots, nt): sum of squares for each partition per day
################################################################
SumSquares <- function(y.by.knots, mu.y, prec, partition, 
                       nt, nknots){

  ss <- matrix(NA, nknots, nt)

  for (k in 1:nknots) {
    these  <- which(partition == k) # identify sites to includes
    ns.k   <- length(these)
    if (ns.k > 0) {
      res.k  <- y.by.knots[[k]] - mu.y[[k]]
      prec.k <- prec[[k]]
      for (t in 1:nt) {
        ss[k, t] <- t(res.k[, t]) %*% prec.k %*% res.k[, t]
      }
    } else {
      ss[k, ] <- rep(0, nt)
    }
  }
    
  return(ss)
}

################################################################
# Arguments:
#   prec.beta(p, p): prior precision of beta
#   e.beta(p): prior mean of beta
#   x.by.knots[nknots]: list of covariate arrays(ns.k, nt, p)
#   y.by.knots[nknots]: list of observed data matrices(ns.k, nt)
#   prec[nknots]: list of precision matrices(ns.k, ns.k)
#   partition(ns): vector of partition membership
#   r.inv(nknots, nt): matrix of random effects at each partition (var scale)
#   nt(1): number of days
#
# Returns:
#   list: 
#     vvv(p, p): posterior variance of beta
#     mmm(p): posterior mean of beta
################################################################
BetaPosterior <- function(prec.beta, e.beta, x.by.knots, y.by.knots, 
						  prec, partition, r.inv, nt){
  
  nknots <- nrow(r.inv)
  p      <- length(e.beta)
  vvv    <- prec.beta
  mmm    <- e.beta

  for (k in 1:nknots) {
    these   <- which(partition == k) # identify sites to includes
    ns.k    <- length(these)
    if (ns.k > 0) {
      x.k     <- x.by.knots[[k]]
      prec.k  <- prec[[k]]
      y.k     <- y.by.knots[[k]] 
      r.inv.k <- r.inv[k, ]
     
      for(t in 1:nt){
        x.kt <- matrix(x.k[, t, ], ns.k, p)
        ttt  <- t(x.kt) %*% prec.k * r.inv.k[t]
        vvv  <- vvv + ttt %*% x.kt
        mmm  <- mmm + ttt %*% y.k[, t]
      }
    }
  }
	
  vvv <- chol2inv(chol(vvv))

  results <- list(vvv=vvv, mmm=mmm)

  return(results)
}

################################################################
# Arguments:
#   ss(nknots, nt): sum of squares for each partition per da
#   log.det(nknots): logdet(prec)
#   r.inv(nknots, nt): matrix of random effects at each partition (var scale)
#   partition(ns): vector of partition membership
#
# Returns:
#   llike(nknots, nt): (log)likelihood
################################################################
LLike <- function(ss, log.det, r.inv, partition, log=TRUE){
  nknots <- nrow(r.inv)
  log.like  <- matrix(NA, nknots, nt)
    
  for(k in 1:nknots){
    these <- which(partition == k) # identify sites to includes
    ns.k  <- length(these)
    r.inv.k <- r.inv[k, ]
    log.like[k,] <- 0.5 * log.det[k] + 0.5 * ns.k * (log(r.inv.k)) - 
                    0.5 * ss[k, ] * r.inv.k
  }
    
  if(!log){
    log.like <- exp(log.like)
  }
    
  return(log.like)
}

################################################################
# Arguments:
#   preds(yp, nt, iters): mcmc predictions at validation
#                         locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np, nt): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
QuantScore <- function(preds, probs, validate){
  nt <- ncol(validate)
  np <- nrow(validate)
  nprobs <- length(probs)
        
  # apply gives nprobs x nsites. looking to find each site's quantile over all
  # of the days.
  pred.quants <- apply(preds, 1, quantile, probs=probs, na.rm=T)
    
  scores.sites <- array(NA, dim=c(nprobs, np, nt))
    
  for (q in 1:nprobs) {
    diff <- pred.quants[q, ] - validate
    i <- ifelse(diff >= 0, 1, 0)
    scores.sites[q, , ] <- 2 * (i - probs[q]) * diff
  }
    
  scores <- apply(scores.sites, 1, mean, na.rm=T)

  return(scores)
}

################################################################
# Arguments:
#   preds(yp, nt, iters): mcmc predictions at validation
#                         locations
#   probs(nthreshs): sample quantiles for scoring
#   validate(np, nt): validation data
#
# Returns:
#   list:
#     scores(nthreshs): a single brier score per threshold
#     threshs(nthreshs): sample quantiles from dataset
################################################################
BrierScore <- function(preds, probs, validate){
  nthreshs <- length(probs)
  thresholds <- quantile(validate, probs=probs, na.rm=T)
    
  scores <- rep(NA, nthreshs)
  for (b in 1:nthreshs) {
    pat <- apply((preds > thresholds[b]), c(1, 2), mean)
    ind <- validate < thresholds[b]
    scores[b] <- mean((ind - pat)^2, na.rm=T)
  }
    
  return(scores)
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