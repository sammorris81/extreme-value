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

dfoldnorm <- function(x, mu, sig) {
  d <- dnorm(x, mu, sig) + dnorm(x, -mu, sig)
  return(d)
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

#### Go from normal to Gamma(alpha, beta)
# Arguments:
#   tau.star(nknots, nt): copula terms for each knot / day
#   phi(1): AR(1) coefficient
#   alpha(1): gamma shape parameter
#   beta(1): gamma rate parameter
cop.IG <- function(tau.star, phi, alpha, beta) {

  if (!is.matrix(tau.star)) {
    stop("Error cop.IG: Must have a matrix for the tau.star terms")
  }

  nt <- ncol(tau.star)
  nknots <- nrow(tau.star)

  mean <- matrix(0, nknots, nt)  # first day is mean 0
  for (t in 2:nt) {
    mean[, t] <- phi * tau.star[, (t-1)]
  }

  sd <- sqrt(1 - phi^2)
  res.std <- (tau.star - mean) / sd
  tau <- qgamma(pnorm(res.std), shape=alpha, rate=beta)

  return(tau)
}

get.tau.mh.idx <- function(nparts, ns, mh.tau.parts) {
  idx <- max(which((nparts / ns) >= mh.tau.parts))
  return(idx)
}

#### Go from Gamma(alpha, beta) to normal
# Arguments:
#   tau(nknots, nt): variance terms for each knot / day
#   phi(1): AR(1) coefficient
#   alpha(1): gamma shape parameter
#   beta(1): gamma rate parameter
cop.inv.IG <- function(tau, phi, alpha, beta) {

  if (!is.matrix(tau)) {
    stop("Error cop.inv.IG: Must have a matrix for the tau terms")
  }

  nt <- ncol(tau)
  nknots <- nrow(tau)

  tau.star <- matrix(0, nknots, nt)
  for (t in 1:nt) {
    if (t == 1) {
      mean <- 0
    } else {
      mean <- phi * tau.star[, (t - 1)]
    }
    sd <- sqrt(1 - phi^2)
    tau.star[, t] <- qnorm(pgamma(tau[, t], shape=alpha, rate=beta), mean=mean, sd=sd)
  }

  return(tau.star)
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


rpotspat <- function(nt, x, s, beta, gamma, nu, gau.rho, t.rho,
                     mixprob, lambda, tau.alpha, tau.beta, nknots) {

  p <- dim(x)[3]
  ns <- nrow(s)

  y <- matrix(NA, ns, nt)
  tau <- matrix(NA, nknots, nt)
  z <- matrix(NA, nknots, nt)
  g <- matrix(NA, ns, nt)

  d <- as.matrix(dist(s))
  # gau is used if mixprob = 0
  gau.C      <- CorFx(d=d, gamma=gamma, rho=gau.rho, nu=nu)
  gau.tau    <- matrix(0.25, nrow=nknots, ncol=nt)
  gau.sd     <- 1 / sqrt(gau.tau)
  gau.z      <- gau.sd * matrix(abs(rnorm(nknots * nt, 0, 1)), nknots, nt)

  # t is used if mixprob = 1
  t.C      <- CorFx(d=d, gamma=gamma, rho=t.rho, nu=nu)
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
      tau[, t] <- t.tau[, t]
      taug     <- t.tau[g, t]
      z[, t]   <- t.z[, t]
      zg       <- t.z[g, t]
      C        <- t.C
    } else {
      tau[, t] <- gau.tau[, t]
      taug     <- gau.tau[g, t]
      z[, t]   <- gau.z[, t]
      zg       <- gau.z[g, t]
      C        <- gau.C
    }

    sdg  <- 1 / sqrt(taug)
    C <- diag(sdg) %*% C %*% diag(sdg)
    chol.C <- chol(C)

    if (p == 1) {
      x.beta <- matrix(x[, t, ], ns, 1) * beta
    } else {
      x.beta <- x[, t, ] %*% beta
    }
    mu <- x.beta + lambda * zg

    y.t <- mu + t(chol.C) %*% matrix(rnorm(ns), ns, 1)
    y[, t] <- y.t
  }

  results <- list(y=y, tau=tau, z=z, knots=knots)
}

makeknotsTS <- function(nt, nknots, s, phi) {
  knots <- array(NA, dim=c(nknots, nt, 2))
  min.s1 <- min(s[, 1])
  max.s1 <- max(s[, 1])
  range.s1 <- max.s1 - min.s1
  min.s2 <- min(s[, 2])
  max.s2 <- max(s[, 2])
  range.s2 <- max.s2 - min.s2

  # starting point is uniform over the space
  knots[, 1, 1] <- runif(nknots, min.s1, max.s1)
  knots[, 1, 2] <- runif(nknots, min.s2, max.s2)

  for (t in 2:nt) {
    # transform the previous day's knots to (0, 1)
    knots.u.lag1 <- (knots[, (t - 1), 1] - min.s1) / range.s1
    knots.star.lag1 <- qnorm(knots.u.lag1)
    knots.u.lag2 <- (knots[, (t - 1), 2] - min.s2) / range.s2
    knots.star.lag2 <- qnorm(knots.u.lag2)

    # draw the new set of knots using the time series
    knots.star.1  <- phi * knots.star.lag1 + sqrt(1 - phi^2) * rnorm(nknots)
    knots.u.1     <- pnorm(knots.star.1)
    knots[, t, 1] <- knots.u.1 * range.s1 + min.s1
    knots.star.2  <- phi * knots.star.lag2 + sqrt(1 - phi^2) * rnorm(nknots)
    knots.u.2     <- pnorm(knots.star.2)
    knots[, t, 2] <- knots.u.2 * range.s2 + min.s2
  }

  return(knots)
}

maketauTS <- function(nt, nknots, tau.alpha, tau.beta, phi) {
  tau.star <- matrix(NA, nrow=nknots, ncol=nt)
  tau.star[, 1] <- rnorm(nknots, 0, 1)
  for (t in 2:nt) {
    tau.star[, t] <- phi * tau.star[, (t - 1)] + sqrt(1 - phi^2) * rnorm(nknots)
  }

  tau <- qgamma(pnorm(tau.star), tau.alpha, tau.beta)
  return(tau)
}

makezTS <- function(nt, nknots, tau, phi) {
  sigma <- 1 / sqrt(tau)
  z.star <- matrix(NA, nrow=nknots, ncol=nt)
  z.star[, 1] <- rnorm(nknots, 0, sqrt(sigma[, 1]))
  for (t in 2:nt) {
    z.star[, t] <- phi * z.star[, (t - 1)] + sqrt(sigma[, t] * (1 - phi^2)) * rnorm(nknots)
  }

  z <- abs(z.star)
  return(z)
}

rpotspatTS <- function(nt, x, s, beta, gamma, nu, rho, phi.z, phi.w, phi.tau,
                       lambda, tau.alpha, tau.beta, nknots) {

  p <- dim(x)[3]
  ns <- nrow(s)

  y <- matrix(NA, ns, nt)
  tau <- matrix(NA, nknots, nt)
  z <- matrix(NA, nknots, nt)
  g <- matrix(NA, ns, nt)

  d <- as.matrix(dist(s))

  C      <- CorFx(d=d, gamma=gamma, rho=rho, nu=nu)
  tau    <- maketauTS(nt=nt, nknots=nknots, tau.alpha=tau.alpha,
                      tau.beta=tau.beta, phi=phi.tau)
  sd     <- 1 / sqrt(tau)
  z      <- makezTS(nt=nt, nknots=nknots, tau=tau, phi=phi.z)

  knots <- makeknotsTS(nt=nt, nknots=nknots, s=s, phi=phi.w)

  for (t in 1:nt) {
    knots.t <- matrix(knots[, t, ], nknots, 2)
    g       <- mem(s, knots.t)
    zg      <- z[g, t]
    taug    <- tau[g, t]

    sdg  <- 1 / sqrt(taug)
    C <- diag(sdg) %*% C %*% diag(sdg)
    chol.C <- chol(C)

    if (p == 1) {
      x.beta <- x[, t, , drop=F] * beta
    } else {
      x.beta <- x[, t, ] %*% beta
    }
    mu <- x.beta + lambda * zg

    y.t <- mu + t(chol.C) %*% matrix(rnorm(ns), ns, 1)
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