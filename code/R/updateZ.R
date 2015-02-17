# we need to include zg for the multiple knot update
updateZ <- function(y, x.beta, zg, prec, tau, mu, taug, g, lambda.1, lambda.2) {
  nknots <- nrow(tau)
  nt     <- ncol(tau)
  ns     <- nrow(y)

  if (nknots == 1) {
    z <- matrix(NA, nknots, nt)
    zg <- matrix(NA, ns, nt)
    for (t in 1:nt) {
      res.t <- y[, t] - x.beta[, t]
      mmm   <- lambda.1 * tau[1, t] * sum(prec %*% res.t)
      vvv   <- tau[1, t] * (lambda.2 + sum(prec))

      vvv <- 1 / vvv
      mmm <- vvv * mmm

      z[1, t] <- abs(rnorm(1, mmm, sqrt(vvv)))
      zg[, t] <- z[1, t]
    }
  } else {
    z.update <- z.Rcpp(taug=taug, tau=tau, y=y, x_beta=x.beta, mu=mu, g=g,
                       prec=prec, lambda_1=lambda.1, lambda_2=lambda.2,
                       zg=zg)
    z  <- z.update$z
    zg <- z.update$zg
  }
  results <- list(z=z, zg=zg)
  return(results)
}

updateZTS <- function(z, zg, y, lambda.1, lambda.2, x.beta,
                      phi, tau, taug, g, prec,
                      acc, att, mh, acc.phi, att.phi, mh.phi) {
  nt <- ncol(z)
  nknots <- nrow(z)

  # transform via copula to normal
  sig    <- 1 / sqrt(tau * lambda.2)
  z.star <- hn.cop(x=z, sig=sig)

  # storage for block update
  # z.new      <- z
  # zg.new     <- zg

  for (t in 1:nt) {
  	taug.t  <- sqrt(taug[, t])
    mu.t    <- x.beta[, t] + lambda.1 * zg[, t]
    cur.res <- y[, t] - mu.t
    cur.lly <- -0.5 * quad.form(prec, taug.t * cur.res)

    for (k in 1:nknots) {
      att[k, t]     <- att[k, t] + 1
      these         <- which(g[, t] == k)
      can.z.star    <- z.star[, t]
      can.z.star[k] <- rnorm(1, z.star[k, t], mh[k, t])

      # transform back to R+
      can.z  <- hn.invcop(x=can.z.star, sig=sig[, t])
      if (can.z[k] < 1e-6) {  # numerical stability
        can.z[k] <- 1e-6
      }
      can.zg <- can.z[g[, t]]  # ns long

      can.mu.t <- x.beta[, t] + lambda.1 * can.zg
      can.res  <- y[, t] - can.mu.t
      can.lly  <- -0.5 * quad.form(prec, taug.t * can.res)

      # prior
      if (t > 1) {
        mean <- phi * z.star[k, (t - 1)]
        sd   <- sqrt(1 - phi^2)
      } else {
        mean <- 0
        sd   <- 1
      }

      R <- can.lly - cur.lly +
           dnorm(can.z.star[k], mean, sd, log=T) -
           dnorm(z.star[k, t], mean, sd, log=T)

      if (t < nt) {
      	sd <- sqrt(1 - phi^2)
        can.mean <- phi * can.z.star[k]
        cur.mean <- phi * z.star[k, t]
        z.star.lag1 <- z.star[k, t + 1]
      	R <- R + dnorm(z.star.lag1, can.mean, sd, log=T) -
                 dnorm(z.star.lag1, cur.mean, sd, log=T)
      }

      if (!is.na(R)) { if (log(runif(1)) < R) {
        acc[k, t]    <- acc[k, t] + 1
        z[k, t]      <- can.z[k]
        zg[these, t] <- can.z[k]
        z.star[k, t] <- can.z.star[k]
        cur.lly      <- can.lly
      } }
    }
  }

  # zg     <- zg.new
  # z      <- z.new

  # phi.z
  phi.update <- updatePhiTS(data=z.star, phi=phi, day.mar=2,
                            att=att.phi, acc=acc.phi, mh=mh.phi)
  phi     <- phi.update$phi
  acc.phi <- phi.update$acc
  att.phi <- phi.update$att

  results <- list(z=z, acc=acc, att=att, zg=zg,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi)

  return(results)
}