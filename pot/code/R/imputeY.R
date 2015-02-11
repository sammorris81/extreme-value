# y is data matrix
# thresh.mtx is matrix of upper thresholds
# taug is site-specific precision terms
# mu is mean
# obs is which observations need to be imputed

imputeY <- function(y, taug, mu, obs, cor, gamma, thresh.mtx=NULL) {
  ns <- nrow(y)
  nt <- ncol(y)

  y.impute <- matrix(y, ns, nt)
  res <- y - mu

  # this will not change for each day - only tau does
  # may need to rethink when gamma is close to 0 or 1
  if (gamma < 1e-6) {               # observations are completely independent
    vvv <- diag(1, ns)
  } else if ((1 - gamma) < 1e-6) {  # observations are completely dependent
    vvv <- cor
  } else {                          # observations are somewhere in between
    vvv <- chol2inv(chol(cor)) / gamma + diag(1 / (1 - gamma), ns)
    vvv <- chol2inv(chol(vvv))
  }
  t.chol.vvv <- t(chol(vvv))

  for (t in 1:nt) {
    taug.t <- sqrt(taug[, t])
    mu.t   <- mu[, t]
    res.t  <- res[, t]
    impute.these <- which(obs[, t])

    if (length(impute.these) > 0) {

      # draw theta for the spatial random effect
      mmm   <- vvv %*% (res[, t] / (1 - gamma))
      theta <- mmm + t.chol.vvv %*% rnorm(ns) / taug.t

      # get the expected value for truncated imputation
      impute.e     <- mu.t[impute.these] + theta[impute.these]
      impute.sd    <- sqrt(1 - gamma) / taug.t[impute.these]

      if (is.null(thresh.mtx)) {
        u.upper <- rep(1, length(impute.these))
      } else {
        thresh.these <- thresh.mtx[impute.these, t]
        u.upper <- pnorm(thresh.these, impute.e, impute.sd)
      }

      u.impute   <- runif(length(impute.these))
      y.impute.t <- impute.e + impute.sd * qnorm(u.impute * u.upper)
      y.impute.t[(u.upper < 1e-6)] <- 0.99999 * thresh.these[(u.upper < 1e-6)]
      # y.impute.t <- ifelse(  # for numerical stability
      #   u.upper < 1e-6,
      #   thresh.mtx.fudge[impute.these, t],
      #   impute.e + impute.sd * qnorm(u.impute * u.upper)
      # )
      y.impute[impute.these, t] <- y.impute.t
    }
  }

  return(y.impute)
}

