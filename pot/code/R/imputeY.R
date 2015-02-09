# y is data matrix
# thresh.mtx is matrix of upper thresholds
# taug is site-specific precision terms
# mu is mean
# obs is which observations need to be imputed

imputeY <- function(y, taug, mu, obs, thresh.mtx=NULL) {
  ns <- nrow(y)
  nt <- ncol(y)

  y.impute <- matrix(y, ns, nt)

  for (t in 1:nt) {
    taug.t <- sqrt(taug[, t])
    mu.t   <- mu[, t]
    res.t  <- res[, t]
    impute.these <- which(obs[, t])

    if (length(impute.these) > 0) {
      # c function to find all conditional means and standard deviations
      impute.cond <- conditional.mean(mn=mu.t, prec=prec.cor, res=res.t,
                                      taug=taug.t, include=impute.these)
      impute.sd    <- impute.cond$cond.sd
      impute.e     <- impute.cond$cond.mn
      thresh.these <- thresh.mtx[impute.these, t]

      if (is.null(thresh.mtx)) {
        u.upper <- rep(1, length(impute.these))
      } else {
        thresh.mtx.fudge <- 0.99999 * thresh.mtx  # for numerical stability
        u.upper <- pnorm(thresh.these, impute.e, impute.sd)

      }

      u.impute   <- runif(length(impute.these))
      y.impute.t <- ifelse(  # for numerical stability
        u.upper < 1e-6,
        thresh.mtx.fudge[impute.these, t],
        impute.e + impute.sd * qnorm(u.impute * u.upper)
      )
      y.impute[impute.these, t] <- y.impute.t
    }
  }

  return(y.impute)
}

