# beta.m is prior mean for beta terms
# beta.s is prior sd for beta terms
# x is design matrix
# y is observations
# taug is site-specific precision terms
# prec is inverse correlation matrix
updateBeta <- function(beta.m, beta.s, x, y, zg, lambda.1, taug, prec) {
  p <- dim(x)[3]
  t <- ncol(y)

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