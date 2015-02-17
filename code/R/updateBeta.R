# beta.m is prior mean for beta terms
# beta.s is prior sd for beta terms
# x is design matrix
# y is observations
# taug is site-specific precision terms
# prec is inverse correlation matrix
updateBeta <- function(beta.m, beta.s, x, y, zg, lambda.1, taug, prec) {
  p <- dim(x)[3]
  nt <- ncol(y)

  if (length(beta.m != p)) {
    mmm <- rep(beta.m, p)
  } else {
    mmm <- beta.m
  }

  vvv <- diag(p) / beta.s^2

  for (t in 1:nt) {
    taug.t <- sqrt(taug[, t])
    x.t    <- x[, t, ] * taug.t
    res.t  <- (y[, t] - lambda.1 * zg[, t]) * taug.t
    ttt    <- t(x.t) %*% prec
    vvv    <- vvv + ttt %*% x.t
    mmm    <- mmm + ttt %*% res.t
  }

  vvv  <- chol2inv(chol(vvv))
  mmm  <- vvv %*% mmm
  beta <- mmm + t(chol(vvv)) %*% rnorm(length(mmm))

  return (beta)
}

# another way to get the regression coefficients
# using the parameterization y = xbeta + lambda |z| where z(s) ~ HN(0, sig(s))
# if skew model, it treats z(s) as another covariate in the model and updates
# lambda as a coefficient
updateBeta1 <- function(beta.m, beta.s, x, y, zg, taug, prec, skew,
                        lambda, lambda.m, lambda.s) {
  p <- dim(x)[3]
  nt <- ncol(y)

  if (length(beta.m != p)) {
    mmm <- rep(beta.m, p)
  } else {
    mmm <- beta.m
  }

  if (skew) {
    mmm <- c(mmm, lambda.m)
    vvv <- diag(p + 1) / beta.s^2
    vvv[p + 1, p + 1] <- 1 / lambda.s^2
  } else {
    vvv <- diag(p) / beta.s^2
  }

  for (t in 1:nt) {
    taug.t <- sqrt(taug[, t])
    if (skew) {
      x.t  <- cbind(x[, t, ], zg[, t]) * taug.t
    } else {
      x.t  <- x[, t, ] * taug.t
    }
    y.t    <- y[, t] * taug.t
    ttt    <- t(x.t) %*% prec
    vvv    <- vvv + ttt %*% x.t
    mmm    <- mmm + ttt %*% y.t
  }

  vvv  <- chol2inv(chol(vvv))
  mmm  <- vvv %*% mmm
  beta <- mmm + t(chol(vvv)) %*% rnorm(length(mmm))

  if (skew) {
    results <- list(beta=beta[1:p], lambda=beta[p + 1])
  } else {
    results <- list(beta=beta, lambda=0)
  }

  return(results)
}