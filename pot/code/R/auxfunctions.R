################################################################################
# Common data transformations
################################################################################
transform <- list(
  logit = function(x, lower=0, upper=1) {
    x <- (x - lower) / (upper - lower)
    log(x / (1 - x))
  },
  inv.logit = function(x, lower=0, upper=1) {
    p <- exp(x) / (1 + exp(x))
    p <- p * (upper - lower) + lower
  },
  probit = function(x, lower=0, upper=1) {
    x <- (x - lower) / (upper - lower)
    qnorm(x)
  },
  inv.probit = function(x, lower=0, upper=1) {
    p <- pnorm(x)
    p <- p * (upper - lower) + lower
  },
  log = function(x) log(x),
  exp = function(x) exp(x),
  copula = function(dens) {
    this.dens <- paste("p", dens, sep="")
    function(x, ...) qnorm(do.call(this.dens, args=list(x, ...)))
  },
  inv.copula = function(dens) {
    this.dens <- paste("q", dens, sep="")
    function(x, ...) do.call(this.dens, args=list(pnorm(x), ...))
  }
)

################################################################################
# Useful densities
################################################################################

#################################################################
# Half normal: Z* = |Z| is HN(0, sig) when Z ~ N(0, sig)
# Note: the location parameter must be 0
#################################################################
dhn <- function(x, sig=1, log=FALSE) {
  if (sum(x < 0) > 0) {
    stop("x must be non-negative")
  }
  density <- dnorm(x, 0, sig) / 2

  if (!log) {
    return(density)
  } else {
    return(log(density))
  }
}

phn <- function(x, sig=1, lower.tail=TRUE) {
  if (sum(x < 0) > 0) {
    stop("x must be non-negative")
  }
  p.val <- 2 * pnorm(x, 0, sig) - 1
  if (!lower.tail) {
    p.val <- 1 - p.val
  }
  return(p.val)
}

qhn <- function(p, sig=1, lower.tail=TRUE) {  # for the half normal density, mu = 0
  if (!lower.tail) {
    p <- 1 - p
  }
  quant <- qnorm((0.5 * p + 0.5), 0, sig)
  return(quant)
}

################################################################################
# Useful copulas
################################################################################
gamma.cop    <- transform$copula(dens="gamma")
gamma.invcop <- transform$inv.copula(dens="gamma")
hn.cop    <- transform$copula(dens="hn")
hn.invcop <- transform$inv.copula(dens="hn")

################################################################################
# MCMC Metropolis SD update
################################################################################

#################################################################
# Description:
#   Function to handle all updates for the MH sds.
#
# Arguments:
#   acc: current number of acceptances
#   att: current number of attempts
#   mh: current mh sd
#
# Returns:
#   mh: updated mh sd
#################################################################
mhupdate <- function(acc, att, mh, nattempts=50, lower=0.8, higher=1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < 0.25) & these.update
  these.high   <- (acc.rate > 0.50) & these.update

  mh[these.low]  <- mh[these.low] * lower
  mh[these.high] <- mh[these.high] * higher

  acc[these.update] <- 0
  att[these.update] <- 0

  results <- list(acc=acc, att=att, mh=mh)
  return(results)
}

loglikeY <- function(y, taug, mu, obs, prec, thresh.mtx) {

  if (is.null(dim(y))) {
    ns <- length(y)
    nt <- 1
  } else {
    ns <- nrow(y)
    nt <- ncol(y)
  }

  y.ll <- ifelse(
    y <= thresh.mtx,
    thresh.mtx,
    y
  )
  y.ll <- matrix(y.ll, ns, nt)
  mu   <- matrix(mu, ns, nt)
  obs  <- matrix(obs, ns, nt)
  ll <- matrix(NA, ns, nt)

  for (t in 1:nt) {
    taug.t <- sqrt(taug[, t])
    mu.t   <- mu[, t]
    res.t  <- y.ll[, t] - mu.t

    cond <- conditional.mean(mn=mu.t, prec=prec, res=res.t,
                                    taug=taug.t)
    cond.sd <- cond$cond.sd
    cond.e  <- cond$cond.mn

    ll[, t] <- obs[, t] * pnorm(y.ll[, t], cond.e, cond.sd, log=TRUE) +
               !obs[, t] * dnorm(y.ll[, t], cond.e, cond.sd, log=TRUE)
  }

  loglike.y <- sum(ll)
}

# Function to return y' %*% prec %*% y
rss <- function(prec, y) {
  if (is.null(dim(y))) {
    nt <- 1
    y  <- matrix(y, ncol=1)
  } else {
    nt <- ncol(y)
  }
  results <- rep(0, nt)
  for (t in 1:nt) {  # benchmarking suggests loop is as fast as apply
    results[t] <- quad.form(prec, y[, t])
  }
  return(results)
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

CorFx <- function(d, gamma, rho, nu) {
  if (rho < 1e-6) {
  	n <- nrow(d)
  	cor <- diag(1, nrow=n)
  } else {
    if (nu == 0.5) {
      cor <- gamma * simple.cov.sp(D=d, sp.type="exponential", sp.par=c(1, rho), error.var=0,
                                   smoothness=nu, finescale.var=0)
    } else {
      cor <- tryCatch(simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho), error.var=0,
                                    smoothness=nu, finescale.var=0),
                      warning=function(e) {
                        cat("rho = ", rho, "\n")
                        cat("nu = ", nu, "\n")
                      })
      cor <- gamma * cor
    }

    diag(cor) <- 1
  }

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
  d <- rdist(s, knots)
  g <- g.Rcpp(d=d)
  return(g$g)
}

get.tau.mh.idx <- function(nparts, ns, mh.tau.parts) {
  idx <- max(which((nparts / ns) >= mh.tau.parts))
  return(idx)
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


rpotspat <- function(nt, x, s, beta, gamma, nu, rho, dist, lambda,
                     tau.alpha, tau.beta, nknots) {

  p  <- dim(x)[3]
  ns <- nrow(s)
  y  <- matrix(NA, ns, nt)
  z  <- matrix(NA, nknots, nt)
  g  <- matrix(NA, ns, nt)
  d  <- as.matrix(dist(s))
  C  <- CorFx(d=d, gamma=gamma, rho=rho, nu=nu)

  if (dist == "t") {
    tau <- matrix(rgamma(nknots * nt, tau.alpha, tau.beta), nknots, nt)
  } else {
    tau <- matrix(0.25, nrow=nknots, ncol=nt)
  }
  sd <- 1 / sqrt(tau)
  z  <- sd * matrix(abs(rnorm(nknots * nt, 0, 1)), nknots, nt)

  knots  <- array(NA, dim=c(nknots, 2, nt))
  min.s1 <- min(s[, 1])
  max.s1 <- max(s[, 1])
  min.s2 <- min(s[, 2])
  max.s2 <- max(s[, 2])

  for (t in 1:nt) {
    knots[, 1, t] <- runif(nknots, min.s1, max.s1)
    knots[, 2, t] <- runif(nknots, min.s2, max.s2)
    knots.t       <- matrix(knots[, , t], nknots, 2)

    g    <- mem(s, knots.t)
    taug <- tau[g, t]
    zg   <- z[g, t]

    sdg    <- 1 / sqrt(taug)
    chol.C <- chol(diag(sdg) %*% C %*% diag(sdg))

    if (p == 1) {
      x.beta <- matrix(x[, t, ], ns, 1) * beta
    } else {
      x.beta <- x[, t, ] %*% beta
    }
    mu <- x.beta + lambda * zg

    y.t    <- mu + t(chol.C) %*% matrix(rnorm(ns), ns, 1)
    y[, t] <- y.t
  }

  results <- list(y=y, tau=tau, z=z, knots=knots)
}

makeKnotsTS <- function(nt, nknots, s, phi) {
  knots.star <- knots <- array(NA, dim=c(nknots, 2, nt))
  min.s1 <- min(s[, 1])
  max.s1 <- max(s[, 1])
  range.s1 <- max.s1 - min.s1
  min.s2 <- min(s[, 2])
  max.s2 <- max(s[, 2])
  range.s2 <- max.s2 - min.s2

  # starting point is uniform over the space
  knots.star[, 1, 1] <- rnorm(nknots)
  knots.star[, 2, 1] <- rnorm(nknots)

  for (t in 2:nt) {
    # draw the new set of knots using the time series
    knots.star[, 1, t] <- phi * knots.star[, 1, (t - 1)] +
                          sqrt(1 - phi^2) * rnorm(nknots)
    knots.star[, 2, t] <- phi * knots.star[, 2, (t - 1)] +
                          sqrt(1 - phi^2) * rnorm(nknots)
  }
  knots[, 1, ] <- transform$inv.probit(knots.star[, 1, ], lower=min.s1,
                                       upper=max.s1)
  knots[, 2, ] <- transform$inv.probit(knots.star[, 2, ], lower=min.s2,
                                       upper=max.s2)

  return(knots)
}

makeTauTS <- function(nt, nknots, tau.alpha, tau.beta, phi) {
  tau.star <- matrix(NA, nrow=nknots, ncol=nt)
  tau.star[, 1] <- rnorm(nknots, 0, 1)
  for (t in 2:nt) {
    tau.star[, t] <- phi * tau.star[, (t - 1)] + sqrt(1 - phi^2) * rnorm(nknots)
  }

  tau <- gamma.invcop(x=tau.star, tau.alpha, tau.beta)
  return(tau)
}

makeZTS <- function(nt, nknots, tau, phi, lambda.2) {
  z.star <- matrix(NA, nrow=nknots, ncol=nt)
  z.star[, 1] <- rnorm(nknots, 0, 1)
  for (t in 2:nt) {
    z.star[, t] <- phi * z.star[, (t - 1)] + sqrt(1 - phi^2) * rnorm(nknots)
  }

  sd <- 1 / sqrt(tau * lambda.2)
  z  <- hn.invcop(x=z.star, sig=sd)

  return(z)
}

rpotspatTS <- function(nt, x, s, beta, gamma, nu, rho, phi.z, phi.w, phi.tau,
                       lambda, tau.alpha, tau.beta, nknots, dist) {

  p <- dim(x)[3]
  ns <- nrow(s)

  y <- matrix(NA, ns, nt)
  tau <- matrix(NA, nknots, nt)
  z <- matrix(NA, nknots, nt)
  g <- matrix(NA, ns, nt)

  d <- as.matrix(dist(s))

  if (lambda == 0) {
    skew     <- FALSE
    lambda.1 <- 0
  } else {
    skew     <- TRUE
    lambda.1 <- sign(lambda)
    lambda.2 <- 1 / lambda^2
  }

  C <- CorFx(d=d, gamma=gamma, rho=rho, nu=nu)
  if (dist == "t") {
    tau <- makeTauTS(nt=nt, nknots=nknots, tau.alpha=tau.alpha,
                     tau.beta=tau.beta, phi=phi.tau)
  } else if (dist == "gaus") {
    tau <- matrix(0.25, nknots, nt)
  }
  sd <- 1 / sqrt(tau)

  if (skew) {
    z <- makeZTS(nt=nt, nknots=nknots, tau=tau, phi=phi.z, lambda.2=lambda.2)
  } else {
    z <- matrix(0, nrow=nknots, ncol=nt)
  }

  knots <- makeKnotsTS(nt=nt, nknots=nknots, s=s, phi=phi.w)

  for (t in 1:nt) {
    knots.t <- matrix(knots[, , t], nknots, 2)
    g       <- mem(s, knots.t)
    zg      <- z[g, t]
    taug    <- tau[g, t]

    sdg    <- 1 / sqrt(taug)
    C.t    <- diag(sdg) %*% C %*% diag(sdg)
    chol.C <- chol(C.t)

    if (p == 1) {
      x.beta <- x[, t, , drop=F] * beta
    } else {
      x.beta <- x[, t, ] %*% beta
    }

    mu <- x.beta + lambda.1 * zg

    y.t    <- mu + t(chol.C) %*% matrix(rnorm(ns), ns, 1)
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
  pred.quants <- apply(preds, 2, quantile, probs=probs, na.rm=T)  # gives nprobs x np x nt

  scores.sites <- array(NA, dim=c(nprobs, np, nt))

  for (q in 1:nprobs) {
    diff <- pred.quants[q, ] - validate
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