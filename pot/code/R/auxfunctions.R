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
rTNorm <- function(mn, sd, lower=-Inf, upper=Inf, fudge=0) {
  lower.u <- pnorm(lower, mn, sd)
  upper.u <- pnorm(upper, mn, sd)
  
  # replace <- ((mn / sd) > 5) & (lower == 0)
  # lower.u[replace] <- 0
  lower.u <- ifelse( mn / sd > 5 & lower == 0, 0, lower.u )
  U <- tryCatch(runif(length(mn), lower.u, upper.u),
                warning=function(e) { 
                  cat("mn =", mn, "\n")
                  cat("sd =", sd, "\n")
                  cat("lower.u =", lower.u, "\n")
                  cat("upper.u =", upper.u, "\n")	
                })
  y <- qnorm(U, mn, sd)
  
  return(y)
}

CorFx <- function(d, alpha, rho, nu) {
  library(geoR)    
  # using cov.spatial instead of matern because it 
  # doesn't use the bessel function unless needed.
  cor       <- alpha * cov.spatial(d, cov.model="matern", cov.pars=c(1, rho), kappa=nu)
  diag(cor) <- 1

  return(cor)
}

eig.inv <- function(Q, inv=T, logdet=T, mtx.sqrt=T, thresh=0.0000001){
  cor.inv <- NULL
  logdet.prec <- NULL
  cor.sqrt <- NULL
  
  eig <- eigen(Q)
  V <- eig$vectors
  D <- ifelse(eig$values < thresh, thresh, eig$values)
  D.inv <- 1 / D
  
  if (logdet) { logdet.prec <- -0.5 * sum(log(D)) }
  if (inv) { cor.inv <- sweep(V, 2, D.inv, "*") %*% t(V) }
  if (mtx.sqrt) { cor.sqrt <- sweep(V, 2, sqrt(D), "*") %*% t(V) }
  
  results <- list(prec=cor.inv, logdet.prec=logdet.prec, sd.mtx=cor.sqrt)
  
  return(results)
}

chol.inv <- function(Q, inv=T, logdet=T) {
  cor.inv <- NULL
  logdet.prec <- NULL
  chol.Q <- chol(Q)
  
  if (inv) { cor.inv <- chol2inv(chol.Q) }
  if (logdet) { logdet.prec <- -sum(log(diag(chol.Q))) }
  
  results <- list(prec=cor.inv, logdet.prec=logdet.prec, sd.mtx=chol.Q)
  return(results)
}

mem <- function(s, knots) {
  library(fields)
  d <- rdist(s, knots)
  g <- apply(d, 1, which.min)

  return(g)
}

#########################################################################
# Arguments:
#   alpha(1): percentage of variation from spatial
#   lambda(n): eigenvalues from correlation matrix
#
# Returns:
#   logdet.prec(1): logdet of the prec(corr) matrix
#########################################################################
logdet.exp <- function(alpha, lambda) {
  logdet.prec <- -sum(log(1 - alpha + alpha * lambda))
  
  return(logdet.prec)
}



rpotspat <- function(nt, x, s, beta, alpha, nu, gau.rho, t.rho,  
                     mixprob, z.alpha, tau.alpha, tau.beta, nknots) {

  p <- dim(x)[3]
  ns <- nrow(s)
  
  y <- matrix(NA, ns, nt)
  tau <- matrix(NA, nknots, nt)
  z <- matrix(NA, nknots, nt)
  g <- matrix(NA, ns, nt)
   
  d <- as.matrix(dist(s))
  # gau is used if mixprob = 0
  gau.C      <- CorFx(d=d, alpha=alpha, rho=gau.rho, nu=nu)
  gau.C.chol <- chol(gau.C)
  gau.tau    <- matrix(0.25, nrow=nknots, ncol=nt)
  gau.sd     <- 1 / sqrt(gau.tau)
  gau.z      <- gau.sd * matrix(abs(rnorm(nknots * nt, 0, 1)), nknots, nt)
  
  # t is used if mixprob = 1
  t.C      <- CorFx(d=d, alpha=alpha, rho=t.rho, nu=nu)
  t.C.chol <- chol(t.C)
  t.tau    <- matrix(rgamma(nknots * nt, tau.alpha, tau.beta), nknots, nt)
  t.sd     <- 1 / sqrt(t.tau)
  t.z      <- t.sd * matrix(abs(rnorm(nknots * nt, 0, 1)), nknots, nt)
  
  knots <- array(NA, dim=c(nknots, nt, 2))
  min.s1 <- min(s[, 1]); max.s1 <- max(s[, 1])
  min.s2 <- min(s[, 2]); max.s2 <- max(s[, 2])

  for (t in 1:nt) {
    knots[, t, 1] <- runif(nknots, min.s1, max.s1)
    knots[, t, 2] <- runif(nknots, min.s2, max.s2)
    knots.t <- matrix(knots[, t, ], nknots, 2)
    g <- mem(s, knots.t)

    dist <- rbinom(1, 1, mixprob)  # 0: gaussian, 1: t
    if (dist) {
      taug   <- t.tau[g, t]
      zg     <- t.z[g, t]
      chol.C <- gau.C.chol
    } else {
      taug   <- gau.tau[g, t]
      zg     <- gau.z[g, t]
      chol.C <- t.C.chol
    }  
        
    if (p == 1) {
      x.beta <- matrix(x[, t, ], ns, 1) * beta 
    } else {
      x.beta <- x[, t, ] %*% beta
    }
    mu <- x.beta + z.alpha * zg
    
    sdg  <- 1 / sqrt(taug)
    y.t <- t(chol.C) %*% matrix(rnorm(ns), ns, 1)
    y.t <- mu + sdg * y.t
    y[, t] <- y.t
  }
  
  results <- list(y=y, tau=tau, z=z, knots=knots)
}

rpotspatTS <- function(nt, x, s, beta, alpha, nu, gau.rho, t.rho, phi.z, phi.w,
                       mixprob, z.alpha, tau.alpha, tau.beta, nknots) {

  p <- dim(x)[3]
  ns <- nrow(s)
  
  y <- matrix(NA, ns, nt)
  tau <- matrix(NA, nknots, nt)
  z <- matrix(NA, nknots, nt)
  g <- matrix(NA, ns, nt)
   
  d <- as.matrix(dist(s))
  # gau is used if mixprob = 0
  gau.C      <- CorFx(d=d, alpha=alpha, rho=gau.rho, nu=nu)
  gau.C.chol <- chol(gau.C)
  gau.tau    <- matrix(0.25, nrow=nknots, ncol=nt)
  gau.sd     <- 1 / sqrt(gau.tau)
  gau.z      <- gau.sd * matrix(abs(rnorm(nknots * nt, 0, 1)), nknots, nt)
  
  # t is used if mixprob = 1
  t.C      <- CorFx(d=d, alpha=alpha, rho=t.rho, nu=nu)
  t.C.chol <- chol(t.C)
  t.tau    <- matrix(rgamma(nknots * nt, tau.alpha, tau.beta), nknots, nt)
  t.sd     <- 1 / sqrt(t.tau)
  t.z      <- t.sd * matrix(abs(rnorm(nknots * nt, 0, 1)), nknots, nt)
  
  knots <- array(NA, dim=c(nknots, nt, 2))
  min.s1 <- min(s[, 1]); max.s1 <- max(s[, 1])
  min.s2 <- min(s[, 2]); max.s2 <- max(s[, 2])

  for (t in 1:nt) {
    knots[, t, 1] <- runif(nknots, min.s1, max.s1)
    knots[, t, 2] <- runif(nknots, min.s2, max.s2)
    knots.t <- matrix(knots[, t, ], nknots, 2)
    g <- mem(s, knots.t)

    dist <- rbinom(1, 1, mixprob)  # 0: gaussian, 1: t
    if (dist) {
      taug   <- t.tau[g, t]
      zg     <- t.z[g, t]
      chol.C <- gau.C.chol
    } else {
      taug   <- gau.tau[g, t]
      zg     <- gau.z[g, t]
      chol.C <- t.C.chol
    }  
        
    if (p == 1) {
      x.beta <- matrix(x[, t, ], ns, 1) * beta 
    } else {
      x.beta <- x[, t, ] %*% beta
    }
    mu <- x.beta + z.alpha * zg
    
    sdg  <- 1 / sqrt(taug)
    y.t <- t(chol.C) %*% matrix(rnorm(ns), ns, 1)
    y.t <- mu + sdg * y.t
    y[, t] <- y.t
  }
  
  results <- list(y=y, tau=tau, z=z, knots=knots)
}

################################################################
# Arguments:
#   preds(iters, yp, nt): mcmc predictions at validation
#                         locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np, nt): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
QuantScore <- function(preds, probs, validate) {
  nt <- ncol(validate)  # number of prediction days
  np <- nrow(validate)  # number of prediction sites
  nprobs <- length(probs)  # number of quantiles to find quantile score
  
  # we need to know the predicted quantiles for each site and day in the validation set
  pred.quants <- apply(preds, c(2, 3), quantile, probs=probs, na.rm=T)  # gives nprobs x np x nt
  
  scores.sites <- array(NA, dim=c(nprobs, np, nt))
  
  for (q in 1:nprobs) {
    diff <- pred.quants[q, , ] - validate
    i <- diff >= 0  # diff >= 0 means qhat is larger
    scores.sites[q, , ] <- 2 * (i - probs[q]) * diff
  }
  
  scores <- apply(scores.sites, 1, mean, na.rm=T)
  
  return(scores)
}

################################################################
# Arguments:
#   preds(iters, yp, nt): mcmc predictions at validation
#                         locations
#   thresholds(nthreshs): sample quantiles for scoring
#   validate(np, nt): validation data
#
# Returns:
#   list:
#     scores(nthreshs): a single brier score per threshold
#     threshs(nthreshs): sample quantiles from dataset
################################################################
BrierScore <- function(preds, thresholds, validate) {
  nthreshs <- length(thresholds)
  
  scores <- rep(NA, nthreshs)
  for (b in 1:nthreshs) {
    pat <- apply((preds > thresholds[b]), c(2, 3), mean)
    i <- validate > thresholds[b]
    scores[b] <- mean((i - pat)^2, na.rm=T)
  }
  
  return(scores)
}


# QuantScore <- function(preds, probs, validate){
  # nt <- ncol(validate)
  # np <- nrow(validate)
  # nprobs <- length(probs)
        
  # # apply gives nprobs x nsites. looking to find each site's quantile over all
  # # of the days.
  # pred.quants <- apply(preds, 2, quantile, probs=probs, na.rm=T)
    
  # scores.sites <- array(NA, dim=c(nprobs, np, nt))
    
  # for (q in 1:nprobs) {
    # diff <- pred.quants[q, ] - validate
    # i <- ifelse(diff >= 0, 1, 0)
    # scores.sites[q, , ] <- 2 * (i - probs[q]) * diff
  # }
    
  # scores <- apply(scores.sites, 1, mean, na.rm=T)

  # return(scores)
# }


# BrierScore <- function(preds, probs, validate){
  # nthreshs <- length(probs)
  # thresholds <- quantile(validate, probs=probs, na.rm=T)
    
  # scores <- rep(NA, nthreshs)
  # for (b in 1:nthreshs) {
    # pat <- apply((preds > thresholds[b]), c(2, 3), mean)
    # ind <- validate < thresholds[b]
    # scores[b] <- mean((ind - pat)^2, na.rm=T)
  # }
    
  # return(scores)
# }