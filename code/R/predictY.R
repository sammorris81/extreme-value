predictY <- function(d11, d12, cov.model, rho, nu, gamma, res, beta, tau, taug,
                     z, prec, lambda.1, s.pred, x.pred, knots) {
  np <- nrow(d11)
  nt <- ncol(res)
  ns <- nrow(res)
  nknots <- nrow(tau)
  yp <- matrix(NA, np, nt)

  if (cov.model == "matern") {
    s.11 <- gamma * simple.cov.sp(D=d11, sp.type="matern", sp.par=c(1, rho),
                                  error.var=0, smoothness=nu,
                                  finescale.var=0)
    s.12 <- gamma * simple.cov.sp(D=d12, sp.type="matern", sp.par=c(1, rho),
                                  error.var=0, smoothness=nu,
                                  finescale.var=0)
  } else {
    s.11 <- gamma * matrix(exp(-d11 / rho), np, np)
    s.12 <- gamma * matrix(exp(-d12 / rho), np, ns)
  }
  diag(s.11) <- 1
  s.12.22.inv <- s.12 %*% prec
  corp <- s.11 - s.12.22.inv %*% t(s.12)
  corp.sd.mtx <- tryCatch(chol(corp),  # only want the cholesky factor
                          error = function(e) {
                            eig.inv(corp, inv=F, logdet=F, mtx.sqrt=T)$sd.mtx
                          })

  for (t in 1:nt) {
    xp.beta  <- x.pred[, t, ] %*% beta
    if (nknots == 1) {
      gp <- rep(1, np)
    } else {
      gp <- mem(s.pred, knots[, , t])  # find the right partition
    }
    zgp    <- z[gp, t]
    siggp  <- 1 / sqrt(tau[gp, t])  # get the partition's standard deviation
    taug.t <- sqrt(taug[, t])
    mup <- xp.beta + lambda.1 * zgp + siggp *
      s.12.22.inv %*% (taug.t * res[, t])

    yp[, t] <- mup + siggp * t(corp.sd.mtx) %*% rnorm(np, 0, 1)
  }

  return(yp)
}

predictY1 <- function(d11, d12, cov.model, rho, nu, gamma, res, beta, tau, taug,
                     z, prec, lambda, s.pred, x.pred, knots) {
  np <- nrow(d11)
  nt <- ncol(res)
  ns <- nrow(res)
  nknots <- nrow(tau)
  yp <- matrix(NA, np, nt)

  if (cov.model == "matern") {
    s.11 <- gamma * simple.cov.sp(D=d11, sp.type="matern", sp.par=c(1, rho),
                                  error.var=0, smoothness=nu,
                                  finescale.var=0)
    s.12 <- gamma * simple.cov.sp(D=d12, sp.type="matern", sp.par=c(1, rho),
                                  error.var=0, smoothness=nu,
                                  finescale.var=0)
  } else {
    s.11 <- gamma * matrix(exp(-d11 / rho), np, np)
    s.12 <- gamma * matrix(exp(-d12 / rho), np, ns)
  }
  diag(s.11) <- 1
  s.12.22.inv <- s.12 %*% prec
  corp <- s.11 - s.12.22.inv %*% t(s.12)
  corp.sd.mtx <- tryCatch(chol(corp),  # only want the cholesky factor
                          error = function(e) {
                            eig.inv(corp, inv=F, logdet=F, mtx.sqrt=T)$sd.mtx
                          })
  t.corp.sd.mtx <- t(corp.sd.mtx)

  for (t in 1:nt) {
    xp.beta  <- x.pred[, t, ] %*% beta
    if (nknots == 1) {
      gp <- rep(1, np)
    } else {
      gp <- mem(s.pred, knots[, , t])  # find the right partition
    }
    zgp    <- z[gp, t]
    siggp  <- 1 / sqrt(tau[gp, t])  # get the partition's standard deviation
    taug.t <- sqrt(taug[, t])
    mup <- xp.beta + lambda * zgp + siggp * s.12.22.inv %*% (taug.t * res[, t])

    yp[, t] <- mup + siggp * t.corp.sd.mtx %*% rnorm(np, 0, 1)
  }

  return(yp)
}