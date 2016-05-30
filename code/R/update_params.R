# we need to include zg for the multiple knot update
updateZ_disc_lambda <- function(y, x.beta, zg, prec, tau, mu, taug, g,
                                lambda.1, lambda.2) {
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
    z.update <- z.Rcpp_disc_lambda(taug = taug, tau = tau, y = y,
                                   x_beta = x.beta, mu = mu, g = g, prec = prec,
                                   lambda_1 = lambda.1, lambda_2 = lambda.2,
                                   zg = zg)
    z  <- z.update$z
    zg <- z.update$zg
  }
  results <- list(z = z, zg = zg)
  return(results)
}


updateZ_cont_lambda <- function(y, x.beta, zg, prec, tau, mu, taug, g, lambda) {
  nknots <- nrow(tau)
  nt     <- ncol(tau)
  ns     <- nrow(y)

  if (nknots == 1) {
    z <- matrix(NA, nknots, nt)
    zg <- matrix(NA, ns, nt)
    for (t in 1:nt) {
      res.t <- y[, t] - x.beta[, t]
      mmm   <- lambda * sum(prec %*% res.t)
      vvv   <- (1 + lambda^2 * sum(prec))

      vvv <- 1 / vvv
      mmm <- vvv * mmm

      z[1, t] <- abs(rnorm(1, mmm, sqrt(vvv / tau[1, t])))
      zg[, t] <- z[1, t]
    }
  } else {
    z.update <- z.Rcpp_cont_lambda(taug = taug, tau = tau, y = y,
                                   x_beta = x.beta, mu = mu, g = g,
                                   prec = prec, lambda = lambda, zg = zg)
    z  <- z.update$z
    zg <- z.update$zg
  }
  results <- list(z = z, zg = zg)
  return(results)
}

updateZTS_disc_lambda <- function(z, zg, y, lambda.1, lambda.2, x.beta,
                                  phi, tau, taug, g, prec,
                                  acc, att, mh, acc.phi, att.phi, mh.phi) {
  nt <- ncol(z)
  nknots <- nrow(z)

  # transform via copula to normal
  sig    <- 1 / sqrt(tau * lambda.2)
  z.star <- hn.cop(x = z, sig = sig)

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
      can.z  <- hn.invcop(x = can.z.star, sig = sig[, t])
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
        dnorm(can.z.star[k], mean, sd, log = TRUE) -
        dnorm(z.star[k, t], mean, sd, log = TRUE)

      if (t < nt) {
        sd <- sqrt(1 - phi^2)
        can.mean <- phi * can.z.star[k]
        cur.mean <- phi * z.star[k, t]
        z.star.lag1 <- z.star[k, t + 1]
        R <- R + dnorm(z.star.lag1, can.mean, sd, log = TRUE) -
          dnorm(z.star.lag1, cur.mean, sd, log = TRUE)
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

  # phi.z
  phi.update <- updatePhiTS(data = z.star, phi = phi, day.mar = 2,
                            att = att.phi, acc = acc.phi, mh = mh.phi)
  phi     <- phi.update$phi
  acc.phi <- phi.update$acc
  att.phi <- phi.update$att

  results <- list(z = z, acc = acc, att = att, zg = zg,
                  phi = phi, att.phi = att.phi, acc.phi = acc.phi)

  return(results)
}

updateZTS_cont_lambda <- function(z, zg, y, lambda, x.beta,
                                 phi, tau, taug, g, prec,
                                 acc, att, mh, acc.phi, att.phi, mh.phi) {
  nt <- ncol(z)
  nknots <- nrow(z)

  # transform via copula to normal
  sig    <- 1 / sqrt(tau)
  z.star <- hn.cop(x = z, sig = sig)

  for (t in 1:nt) {
    taug.t  <- sqrt(taug[, t])
    mu.t    <- x.beta[, t] + lambda * zg[, t]
    cur.res <- y[, t] - mu.t
    cur.lly <- -0.5 * quad.form(prec, taug.t * cur.res)

    for (k in 1:nknots) {
      att[k, t]     <- att[k, t] + 1
      these         <- which(g[, t] == k)
      can.z.star    <- z.star[, t]
      can.z.star[k] <- rnorm(1, z.star[k, t], mh[k, t])

      # transform back to R+
      can.z  <- hn.invcop(x = can.z.star, sig = sig[, t])
      if (can.z[k] < 1e-6) {  # numerical stability
        can.z[k] <- 1e-6
      }
      can.zg <- can.z[g[, t]]  # ns long

      can.mu.t <- x.beta[, t] + lambda * can.zg
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
        dnorm(can.z.star[k], mean, sd, log = TRUE) -
        dnorm(z.star[k, t], mean, sd, log = TRUE)

      if (t < nt) {
        sd <- sqrt(1 - phi^2)
        can.mean <- phi * can.z.star[k]
        cur.mean <- phi * z.star[k, t]
        z.star.lag1 <- z.star[k, t + 1]
        R <- R + dnorm(z.star.lag1, can.mean, sd, log = TRUE) -
          dnorm(z.star.lag1, cur.mean, sd, log = TRUE)
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

  # phi.z
  phi.update <- updatePhiTS(data = z.star, phi = phi, day.mar = 2,
                            att = att.phi, acc = acc.phi, mh = mh.phi)
  phi     <- phi.update$phi
  acc.phi <- phi.update$acc
  att.phi <- phi.update$att

  results <- list(z = z, acc = acc, att = att, zg = zg,
                  phi = phi, att.phi = att.phi, acc.phi = acc.phi)

  return(results)
}


updateTauGaus <- function(res, prec, tau.alpha, tau.beta) {
  ns <- nrow(res)
  nt <- ncol(res)
  this.rss <- sum(rss(prec = prec, y = res))

  aaa  <- tau.alpha + 0.5 * ns * nt
  bbb  <- tau.beta + 0.5 * this.rss
  tau  <- rgamma(1, aaa, bbb)

  return(tau)
}

updateTau_disc_lambda <- function(tau, taug, g, res, nparts.tau, prec, z,
                                  lambda.2, tau.alpha, tau.beta, skew,
                                  att, acc, mh) {
  ns <- nrow(res)
  nt <- ncol(res)
  nknots <- nrow(tau)

  if (nknots == 1) {
    for (t in 1:nt) {
      rss.t <- quad.form(prec, res[, t])

      aaa <- tau.alpha + 0.5 * ns
      bbb <- tau.beta + 0.5 * rss.t
      if (skew) {  # tau is also in the z likelihood
        aaa <- aaa + 0.5
        bbb <- bbb + 0.5 * z[1, t]^2 * lambda.2
      }

      tau[1, t] <- rgamma(1, aaa, bbb)
      taug[, t] <- tau[1, t]
    }
  } else {  # nknots > 1
    for (t in 1:nt) {
      cur.lly <- 0.5 * sum(log(taug[, t])) -
        0.5 * quad.form(prec, sqrt(taug[, t]) * res[, t])

      for (k in 1:nknots) {
        these  <- which(g[, t] == k)
        nparts <- length(these)
        nparts.tau[k, t] <- nparts

        if (nparts == 0) {  # tau update is from prior
          aaa <- tau.alpha
          bbb <- tau.beta
          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2 * lambda.2
          }

          tau[k, t] <- rgamma(1, aaa, bbb)
          if (tau[k, t] < 1e-6) {
            tau[k, t] <- 1e-6
          }
        } else if (nparts == ns) {
          aaa <- tau.alpha + 0.5 * nparts
          bbb <- tau.beta + 0.5 * quad.form(prec, res[, t])

          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2 * lambda.2
          }

          tau[k, t] <- rgamma(1, aaa, bbb)
          taug[, t] <- tau[k, t]

        } else {  # nparts > 0
          taug.t    <- sqrt(taug[, t])
          att[k, t] <- att[k, t] + 1

          aaa <- 0.5 * nparts
          bbb <- 0.5 * quad.form(prec[these, these], res[these, t])

          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2 * lambda.2
          }

          # the posterior is conjugate when nparts = ns
          # as nparts -> 1, want a wider candidate
          # wider candidate means aaa and bbb decrease
          aaa <- tau.alpha + aaa * nparts / ns
          bbb <- tau.beta + bbb * nparts / ns
          can.tau    <- tau[, t]
          can.tau[k] <- rgamma(1, aaa, bbb)
          if (can.tau[k] < 1e-6) {
            can.tau[k] <- 1e-6
          }
          can.taug   <- can.tau[g[, t]]

          can.lly <- 0.5 * sum(log(can.taug)) -
            0.5 * quad.form(prec, sqrt(can.taug) * res[, t])

          if (skew) {
            cur.llz <- 0.5 * log(tau[k, t]) -
              0.5 * tau[k, t] * z[k, t]^2 * lambda.2
            can.llz <- 0.5 * log(can.tau[k]) -
              0.5 * can.tau[k] * z[k, t]^2 * lambda.2
          } else {
            cur.llz <- can.llz <- 0
          }

          R <- can.lly - cur.lly + can.llz - cur.llz +
            dgamma(can.tau[k], tau.alpha, tau.beta, log = TRUE) -
            dgamma(tau[k, t], tau.alpha, tau.beta, log = TRUE) +
            dgamma(tau[k, t], aaa, bbb, log = TRUE) -
            dgamma(can.tau[k], aaa, bbb, log = TRUE)

          if (!is.na(R)) { if (log(runif(1)) < R) {
            acc[k, t]      <- acc[k, t] + 1
            tau[k, t]      <- can.tau[k]
            taug[these, t] <- can.tau[k]
            cur.lly        <- can.lly
          }}
        }  # fi nparts
      }  # end k
    }  # end t
  }  # fi nknots > 1

  results <- list(tau = tau, taug = taug, acc = acc, att = att)

}

updateTau_cont_lambda <- function(tau, taug, g, res, nparts.tau, prec, z,
                                  tau.alpha, tau.beta, skew, att, acc, mh) {
  ns <- nrow(res)
  nt <- ncol(res)
  nknots <- nrow(tau)

  if (nknots == 1) {
    for (t in 1:nt) {
      rss.t <- quad.form(prec, res[, t])

      aaa <- tau.alpha + 0.5 * ns
      bbb <- tau.beta + 0.5 * rss.t
      if (skew) {  # tau is also in the z likelihood
        aaa <- aaa + 0.5
        bbb <- bbb + 0.5 * z[1, t]^2
      }

      tau[1, t] <- rgamma(1, aaa, bbb)
      taug[, t] <- tau[1, t]
    }
  } else {  # nknots > 1
    for (t in 1:nt) {
      cur.lly <- 0.5 * sum(log(taug[, t])) -
        0.5 * quad.form(prec, sqrt(taug[, t]) * res[, t])

      for (k in 1:nknots) {
        these  <- which(g[, t] == k)
        nparts <- length(these)
        nparts.tau[k, t] <- nparts

        if (nparts == 0) {  # tau update is from prior
          aaa <- tau.alpha
          bbb <- tau.beta
          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2
          }

          tau[k, t] <- rgamma(1, aaa, bbb)
          if (tau[k, t] < 1e-6) {
            tau[k, t] <- 1e-6
          }
        } else if (nparts == ns) {
          aaa <- tau.alpha + 0.5 * nparts
          bbb <- tau.beta + 0.5 * quad.form(prec, res[, t])

          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2
          }

          tau[k, t] <- rgamma(1, aaa, bbb)
          taug[, t] <- tau[k, t]

        } else {  # 0 < nparts < ns
          taug.t    <- sqrt(taug[, t])
          att[k, t] <- att[k, t] + 1

          aaa <- 0.5 * nparts
          bbb <- 0.5 * quad.form(prec[these, these], res[these, t])

          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2
          }

          # the posterior is conjugate when nparts = ns
          # as nparts -> 1, want a wider candidate
          # wider candidate means aaa and bbb decrease
          aaa <- tau.alpha + aaa * nparts / ns
          bbb <- tau.beta + bbb * nparts / ns
          can.tau    <- tau[, t]
          can.tau[k] <- rgamma(1, aaa, bbb)
          if (can.tau[k] < 1e-6) {
            can.tau[k] <- 1e-6
          }
          can.taug   <- can.tau[g[, t]]

          can.lly <- 0.5 * sum(log(can.taug)) -
            0.5 * quad.form(prec, sqrt(can.taug) * res[, t])

          if (skew) {
            cur.llz <- 0.5 * log(tau[k, t]) - 0.5 * tau[k, t] * z[k, t]^2
            can.llz <- 0.5 * log(can.tau[k]) - 0.5 * can.tau[k] * z[k, t]^2
          } else {
            cur.llz <- can.llz <- 0
          }

          R <- can.lly - cur.lly + can.llz - cur.llz +
            dgamma(can.tau[k], tau.alpha, tau.beta, log = TRUE) -
            dgamma(tau[k, t], tau.alpha, tau.beta, log = TRUE) +
            dgamma(tau[k, t], aaa, bbb, log = TRUE) -
            dgamma(can.tau[k], aaa, bbb, log = TRUE)


          if (!is.na(R)) { if (log(runif(1)) < R) {
            acc[k, t]      <- acc[k, t] + 1
            tau[k, t]      <- can.tau[k]
            taug[these, t] <- can.tau[k]
            cur.lly        <- can.lly
          }}
        }  # fi nparts
      }  # end k
    }  # end t
  }  # fi nknots > 1

  results <- list(tau = tau, taug = taug, acc = acc, att = att)

}

updateTauTS_disc_lambda <- function(phi, tau, taug, g, res, nparts.tau, prec,
                                    z, lambda.2, tau.alpha, tau.beta, skew,
                                    att, acc, mh, att.phi, acc.phi, mh.phi) {
  ns <- nrow(res)
  nt <- ncol(res)
  nknots <- nrow(tau)

  tau.star <- gamma.cop(tau, tau.alpha, tau.beta)

  if (nknots == 1) {  # seems to come back ok
    for (t in 1:nt) {
      att[1, t] <- att[1, t] + 1
      cur.lly   <- 0.5 * ns * log(tau[1, t]) -
        0.5 * tau[1, t] * quad.form(prec, res[, t])

      can.tau.star  <- rnorm(1, tau.star[1, t], mh[1, t])

      # transform back to R+
      can.tau <- gamma.invcop(can.tau.star, tau.alpha, tau.beta)
      if (can.tau < 1e-6) {  # numerical stability
        can.tau <- 1e-6
      }

      can.lly <- 0.5 * ns * log(can.tau) -
        0.5 * can.tau * quad.form(prec, res[, t])

      if (skew) {
        cur.llz <- 0.5 * log(tau[1, t]) -
          0.5 * tau[1, t] * z[1, t]^2 * lambda.2
        can.llz <- 0.5 * log(can.tau) -
          0.5 * can.tau * z[1, t]^2 * lambda.2
      } else {
        cur.llz <- can.llz <- 0
      }

      R <- can.lly - cur.lly + can.llz - cur.llz

      if (t > 1) {
        mean <- phi * tau.star[1, (t - 1)]
        sd   <- sqrt(1 - phi^2)
      } else {
        mean <- 0
        sd   <- 1
      }

      # evaluate the prior
      R <- R + dnorm(can.tau.star, mean, sd, log = TRUE) -
        dnorm(tau.star[1, t], mean, sd, log = TRUE)

      # adjust R to account for next day
      if (t < nt) {
        sd            <- sqrt(1 - phi^2)
        can.mean      <- phi * can.tau.star
        cur.mean      <- phi * tau.star[1, t]
        tau.star.lag1 <- tau.star[1, t + 1]

        R <- R + dnorm(tau.star.lag1, can.mean, sd, log = T) -
          dnorm(tau.star.lag1, cur.mean, sd, log = T)
      }

      if (!is.na(R)) { if (log(runif(1)) < R) {
        tau.star[1, t] <- can.tau.star
        tau[1, t]      <- can.tau
        taug[, t]      <- rep(can.tau, ns)
        acc[1, t]      <- acc[1, t] + 1
      }}
    }  # end t
  } else {  # nknots > 1
    for (t in 1:nt) {
      cur.lly <- 0.5 * sum(log(taug[, t])) -
        0.5 * quad.form(prec, sqrt(taug[, t]) * res[, t])

      for (k in 1:nknots) {
        att[k, t] <- att[k, t] + 1
        these  <- which(g[, t] == k)
        nparts <- length(these)
        nparts.tau[k, t] <- nparts

        can.tau.star <- tau.star[, t]
        can.tau.star[k] <- rnorm(1, tau.star[k, t], mh[k, t])

        # transform back to R+
        can.tau  <- gamma.invcop(can.tau.star, tau.alpha, tau.beta)
        if (can.tau[k] < 1e-6) {  # numerical stability
          can.tau[k] <- 1e-6
        }
        can.taug <- can.tau[g[, t]]

        can.lly <- 0.5 * sum(log(can.taug)) -
          0.5 * quad.form(prec, sqrt(can.taug) * res[, t])

        if (skew) {
          cur.llz <- 0.5 * log(tau[k, t]) -
            0.5 * tau[k, t] * z[k, t]^2 * lambda.2
          can.llz <- 0.5 * log(can.tau[k]) -
            0.5 * can.tau[k] * z[k, t]^2 * lambda.2
        } else {
          cur.llz <- can.llz <- 0
        }

        R <- can.lly - cur.lly + can.llz - cur.llz

        # evaluate the prior distribution
        if (t > 1) {
          mean <- phi * tau.star[k, (t - 1)]
          sd   <- sqrt(1 - phi^2)
        } else {
          mean <- 0
          sd   <- 1
        }

        R <- R + dnorm(can.tau.star[k], mean, sd, log = TRUE) -
          dnorm(tau.star[k, t], mean, sd, log = TRUE)

        if (t < nt) {
          sd            <- sqrt(1 - phi^2)
          can.mean      <- phi * can.tau.star[k]
          cur.mean      <- phi * tau.star[k, t]
          tau.star.lag1 <- tau.star[k, t + 1]
          R <- R + dnorm(tau.star.lag1, can.mean, sd, log = TRUE) -
            dnorm(tau.star.lag1, cur.mean, sd, log = TRUE)
        }

        if (!is.na(R)) { if (log(runif(1)) < R) {
          acc[k, t]      <- acc[k, t] + 1
          tau.star[k, t] <- can.tau.star[k]
          tau[k, t]      <- can.tau[k]
          taug[these, t] <- can.tau[k]
          cur.lly        <- can.lly
        }}
      }  # end k
    }  # end t
  }  # fi nknots > 1

  phi.update <- updatePhiTS(data = tau.star, phi = phi, day.mar = 2,
                            att = att.phi, acc = acc.phi, mh = mh.phi)
  phi     <- phi.update$phi
  acc.phi <- phi.update$acc
  att.phi <- phi.update$att

  results <- list(tau = tau, taug = taug, phi = phi, acc = acc, att = att,
                  acc.phi = acc.phi, att.phi = att.phi)
}

updateTauTS_cont_lambda <- function(phi, tau, taug, g, res, nparts.tau, prec,
                                    z, tau.alpha, tau.beta, skew, att, acc, mh,
                                    att.phi, acc.phi, mh.phi) {
  ns <- nrow(res)
  nt <- ncol(res)
  nknots <- nrow(tau)

  tau.star <- gamma.cop(tau, tau.alpha, tau.beta)

  if (nknots == 1) {  # seems to come back ok
    for (t in 1:nt) {
      att[1, t] <- att[1, t] + 1
      cur.lly   <- 0.5 * ns * log(tau[1, t]) -
        0.5 * tau[1, t] * quad.form(prec, res[, t])

      can.tau.star  <- rnorm(1, tau.star[1, t], mh[1, t])

      # transform back to R+
      can.tau <- gamma.invcop(can.tau.star, tau.alpha, tau.beta)
      if (can.tau < 1e-6) {  # numerical stability
        can.tau <- 1e-6
      }

      can.lly <- 0.5 * ns * log(can.tau) -
        0.5 * can.tau * quad.form(prec, res[, t])

      if (skew) {
        cur.llz <- 0.5 * log(tau[1, t]) -
          0.5 * tau[1, t] * z[1, t]^2
        can.llz <- 0.5 * log(can.tau) -
          0.5 * can.tau * z[1, t]^2
      } else {
        cur.llz <- can.llz <- 0
      }

      R <- can.lly - cur.lly + can.llz - cur.llz

      if (t > 1) {
        mean <- phi * tau.star[1, (t - 1)]
        sd   <- sqrt(1 - phi^2)
      } else {
        mean <- 0
        sd   <- 1
      }

      # evaluate the prior
      R <- R + dnorm(can.tau.star, mean, sd, log = TRUE) -
        dnorm(tau.star[1, t], mean, sd, log = TRUE)

      # adjust R to account for next day
      if (t < nt) {
        sd            <- sqrt(1 - phi^2)
        can.mean      <- phi * can.tau.star
        cur.mean      <- phi * tau.star[1, t]
        tau.star.lag1 <- tau.star[1, t + 1]

        R <- R + dnorm(tau.star.lag1, can.mean, sd, log = TRUE) -
          dnorm(tau.star.lag1, cur.mean, sd, log = TRUE)
      }

      if (!is.na(R)) { if (log(runif(1)) < R) {
        tau.star[1, t] <- can.tau.star
        tau[1, t]      <- can.tau
        taug[, t]      <- rep(can.tau, ns)
        acc[1, t]      <- acc[1, t] + 1
      }}
    }  # end t
  } else {  # nknots > 1
    for (t in 1:nt) {
      cur.lly <- 0.5 * sum(log(taug[, t])) -
        0.5 * quad.form(prec, sqrt(taug[, t]) * res[, t])

      for (k in 1:nknots) {
        att[k, t] <- att[k, t] + 1
        these  <- which(g[, t] == k)
        nparts <- length(these)
        nparts.tau[k, t] <- nparts

        can.tau.star <- tau.star[, t]
        can.tau.star[k] <- rnorm(1, tau.star[k, t], mh[k, t])

        # transform back to R+
        can.tau  <- gamma.invcop(can.tau.star, tau.alpha, tau.beta)
        if (can.tau[k] < 1e-6) {  # numerical stability
          can.tau[k] <- 1e-6
        }
        can.taug <- can.tau[g[, t]]

        can.lly <- 0.5 * sum(log(can.taug)) -
          0.5 * quad.form(prec, sqrt(can.taug) * res[, t])

        if (skew) {
          cur.llz <- 0.5 * log(tau[k, t]) -
            0.5 * tau[k, t] * z[k, t]^2
          can.llz <- 0.5 * log(can.tau[k]) -
            0.5 * can.tau[k] * z[k, t]^2
        } else {
          cur.llz <- can.llz <- 0
        }

        R <- can.lly - cur.lly + can.llz - cur.llz

        # evaluate the prior distribution
        if (t > 1) {
          mean <- phi * tau.star[k, (t - 1)]
          sd   <- sqrt(1 - phi^2)
        } else {
          mean <- 0
          sd   <- 1
        }

        R <- R + dnorm(can.tau.star[k], mean, sd, log = TRUE) -
          dnorm(tau.star[k, t], mean, sd, log = TRUE)

        if (t < nt) {
          sd            <- sqrt(1 - phi^2)
          can.mean      <- phi * can.tau.star[k]
          cur.mean      <- phi * tau.star[k, t]
          tau.star.lag1 <- tau.star[k, t + 1]
          R <- R + dnorm(tau.star.lag1, can.mean, sd, log = TRUE) -
            dnorm(tau.star.lag1, cur.mean, sd, log = TRUE)
        }

        if (!is.na(R)) { if (log(runif(1)) < R) {
          acc[k, t]      <- acc[k, t] + 1
          tau.star[k, t] <- can.tau.star[k]
          tau[k, t]      <- can.tau[k]
          taug[these, t] <- can.tau[k]
          cur.lly        <- can.lly
        }}
      }  # end k
    }  # end t
  }  # fi nknots > 1

  phi.update <- updatePhiTS(data = tau.star, phi = phi, day.mar = 2,
                            att = att.phi, acc = acc.phi, mh = mh.phi)
  phi     <- phi.update$phi
  acc.phi <- phi.update$acc
  att.phi <- phi.update$att

  results <- list(tau = tau, taug = taug, phi = phi, acc = acc, att = att,
                  acc.phi = acc.phi, att.phi = att.phi)
}

updateTauAlpha <- function(tau, tau.beta,
                           tau.alpha.min = 0.1, tau.alpha.max = 10,
                           tau.alpha.by = 0.1) {
  # using unique so we can use the same function for all possible setups
  tau <- unique(tau)
  lll <- mmm <- seq(tau.alpha.min, tau.alpha.max, tau.alpha.by)
  for (l in 1:length(lll)) {
    lll[l] <- sum(dgamma(tau, mmm[l], tau.beta, log = TRUE))
  }

  tau.alpha <- sample(mmm, 1, prob = exp(lll - max(lll)))
  return(tau.alpha)
}

updateTauBeta <- function(tau, tau.alpha, tau.beta.a, tau.beta.b) {
  # using unique so we can use the same function for all possible setups
  tau <- unique(tau)
  ntaus <- length(tau)  # typically nknots * nt, but only 1 for gaussian
  a.star <- tau.beta.a + tau.alpha * ntaus
  b.star <- tau.beta.b + sum(tau)

  tau.beta <- rgamma(1, a.star, b.star)
  return(tau.beta)
}

# beta.m is prior mean for beta terms
# beta.s is prior sd for beta terms
# x is design matrix
# y is observations
# taug is site-specific precision terms
# prec is inverse correlation matrix
updateBeta_disc_lambda <- function(beta.m, beta.s, x, y, zg, lambda.1, taug,
                                   prec) {
  p <- dim(x)[3]
  nt <- ncol(y)

  if (length(beta.m != p)) {
    mmm <- rep(beta.m, p)
  } else {
    mmm <- beta.m
  }

  vvv <- diag(p) / (beta.s^2)

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
updateBeta_cont_lambda <- function(beta.m, beta.s, x, y, zg, taug, prec, skew,
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
    vvv <- diag(p + 1) / (beta.s^2)
    vvv[p + 1, p + 1] <- 1 / (lambda.s^2)
  } else {
    vvv <- diag(p) / (beta.s^2)
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
    results <- list(beta = beta[1:p], lambda = beta[p + 1])
  } else {
    results <- list(beta = beta, lambda = 0)
  }

  return(results)
}

updateRhoNu <- function(rho, logrho.m, logrho.s, fixnu, nu, lognu.m, lognu.s,
                        d, gamma, res, taug, prec, cor, logdet.prec, cur.rss,
                        rho.upper = Inf, nu.upper = Inf, att.rho, acc.rho,
                        mh.rho, att.nu, acc.nu, mh.nu) {
  nt <- ncol(res)
  att.rho <- att.rho + 1
  att.nu  <- att.nu + 1

  rho.star     <- transform$probit(rho, lower = 0, upper = rho.upper)
  can.rho.star <- rnorm(1, rho.star, mh.rho)
  can.rho      <- transform$inv.probit(can.rho.star, lower = 0,
                                       upper = rho.upper)

  lognu  <- log(nu)
  if (!fixnu) {
    if (nu.upper == Inf) {
      upper.lognu <- 1
    } else {
      upper.lognu <- pnorm(log(nu.upper), lognu, mh.nu)
    }
    can.nu.u  <- runif(1) * upper.lognu
    can.lognu <- lognu + mh.nu * qnorm(can.nu.u)
  } else {
    can.lognu <- lognu
  }
  can.nu <- exp(can.lognu)

  # CorFx gives can.C too, but want a way to pass through gammaless
  # correlation matrix to speed up gamma update
  can.cor     <- simple.cov.sp(D = d, sp.type = "matern",
                               sp.par = c(1, can.rho),
                               error.var = 0, smoothness = can.nu,
                               finescale.var = 0)
  can.C       <- gamma * can.cor
  diag(can.C) <- 1
  # can.C <- CorFx(d=d, gamma=gamma, rho=can.rho, nu=can.nu)
  can.CC <- tryCatch(chol.inv(can.C, inv = T, logdet = T),
                     error = function(e) {
                       tryCatch(eig.inv(can.C, inv = TRUE, logdet = TRUE,
                                        mtx.sqrt = TRUE),
                                error = function(e) {
                                  print(paste("can.rho =", can.rho))
                                  print(paste("can.nu =", can.nu))
                                })
                     })
  can.prec        <- can.CC$prec
  can.logdet.prec <- can.CC$logdet.prec  # this is the sqrt of logdet.prec

  can.rss <- sum(rss(prec = can.prec, y = sqrt(taug) * res))

  R <- -0.5 * (can.rss - cur.rss) +
    nt * (can.logdet.prec - logdet.prec) +
    dnorm(can.rho.star, log = TRUE) - dnorm(rho.star, log = TRUE)

  if (!fixnu) {
    R <- R + dnorm(can.lognu, lognu.m, lognu.s, log = TRUE) -
      dnorm(lognu, lognu.m, lognu.s, log = TRUE)
    if (upper.lognu < 1) {  # candidate is not symmetric
      R <- R + dnorm(lognu, can.lognu, mh.nu, log = TRUE) -
        pnorm(upper.lognu, can.lognu, mh.nu, log.p = TRUE) -
        dnorm(can.lognu, lognu, mh.nu, log = TRUE) +
        pnorm(upper.lognu, lognu, mh.nu, log.p = TRUE)
    }
  }

  if (!is.na(R)) { if (log(runif(1)) < R) {
    rho         <- can.rho
    nu          <- can.nu
    prec        <- can.prec
    cor         <- can.cor
    logdet.prec <- can.logdet.prec
    acc.rho     <- acc.rho + 1
    acc.nu      <- acc.nu + 1
    cur.rss     <- can.rss
  }}

  results <- list(rho = rho, nu = nu, prec = prec, cor = cor,
                  logdet.prec = logdet.prec, cur.rss = cur.rss,
                  att.rho = att.rho, acc.rho = acc.rho,
                  att.nu = att.nu, acc.nu = acc.nu)
}


updateGamma <- function(gamma, gamma.m, gamma.s, d, rho, nu, taug, res, prec,
                        cor, logdet.prec, cur.rss, att, acc, mh) {
  nt <- ncol(res)
  att <- att + 1

  gamma.star     <- transform$probit(gamma)
  can.gamma.star <- rnorm(1, gamma.star, mh)
  can.gamma      <- transform$inv.probit(can.gamma.star)

  can.C  <- can.gamma * cor
  diag(can.C) <- 1
  can.CC <- tryCatch(chol.inv(can.C, inv = TRUE, logdet = TRUE),
                     error = function(e) {
                       tryCatch(eig.inv(can.C, inv = TRUE, logdet = TRUE,
                                        mtx.sqrt = TRUE),
                                error = function(e) {
                                  print(paste("can.gamma =", can.gamma))
                                })
                     })
  can.prec        <- can.CC$prec
  can.logdet.prec <- can.CC$logdet.prec  # sqrt of logdet.prec

  can.rss <- sum(rss(prec = can.prec, y = sqrt(taug) * res))

  R <- -0.5 * (can.rss - cur.rss) + nt * (can.logdet.prec - logdet.prec) +
    dnorm(can.gamma.star, log = TRUE) - dnorm(gamma.star, log = TRUE)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc         <- acc + 1
    gamma       <- can.gamma
    prec        <- can.prec
    logdet.prec <- can.logdet.prec
    cur.rss     <- can.rss
  }}

  results <- list(gamma = gamma, prec = prec, logdet.prec = logdet.prec,
                  cur.rss = cur.rss, acc = acc, att = att)

  return(results)
}

updateRhoNuGamma <- function(rho, logrho.m, logrho.s, fixnu, nu, lognu.m,
                             lognu.s, d, gamma, res, taug, prec,
                             logdet.prec, cur.rss, rho.upper = Inf,
                             nu.upper = Inf, att.rho, acc.rho, mh.rho,
                             att.nu, acc.nu, mh.nu,
                             att.gamma, acc.gamma, mh.gamma) {
  nt <- ncol(res)
  att.rho <- att.rho + 1
  att.nu  <- att.nu + 1
  att.gamma <- att.gamma + 1

  logrho <- log(rho)
  if (rho.upper == Inf) {
    upper.logrho <- 1
  } else {
    upper.logrho <- pnorm(log(rho.upper), logrho, mh.rho)
  }
  can.rho.u  <- runif(1) * upper.logrho
  can.logrho <- logrho + mh.rho * qnorm(can.rho.u)
  can.rho    <- exp(can.logrho)

  lognu  <- log(nu)
  if (!fixnu) {
    if (nu.upper == Inf) {
      upper.lognu <- 1
    } else {
      upper.lognu <- pnorm(log(nu.upper), lognu, mh.nu)
    }
    can.nu.u  <- runif(1) * upper.lognu
    can.lognu <- lognu + mh.nu * qnorm(can.nu.u)
  } else {
    can.lognu <- lognu
  }
  can.nu <- exp(can.lognu)

  gamma.star     <- transform$probit(gamma)
  can.gamma.star <- rnorm(1, gamma.star, mh.gamma)
  can.gamma      <- transform$inv.probit(can.gamma.star)

  # CorFx gives can.C too, but want a way to pass through gammaless
  # correlation matrix to speed up gamma update
  can.C <- CorFx(d = d, gamma = can.gamma, rho = can.rho, nu = can.nu)
  can.CC <- tryCatch(chol.inv(can.C, inv = TRUE, logdet = TRUE),
                     error = function(e) {
                       tryCatch(eig.inv(can.C, inv = TRUE, logdet = TRUE,
                                        mtx.sqrt = TRUE),
                                error = function(e) {
                                  print(paste("can.rho =", can.rho))
                                  print(paste("can.nu =", can.nu))
                                })
                     })
  can.prec        <- can.CC$prec
  can.logdet.prec <- can.CC$logdet.prec  # this is the sqrt of logdet.prec

  can.rss <- sum(rss(prec = can.prec, y = sqrt(taug) * res))

  R <- -0.5 * (can.rss - cur.rss) +
    nt * (can.logdet.prec - logdet.prec) +
    dnorm(can.logrho, logrho.m, logrho.s, log = TRUE) -
    dnorm(logrho, logrho.m, logrho.s, log = TRUE) +
    dnorm(can.gamma.star, log = TRUE) - dnorm(gamma.star, log = TRUE)

  if (upper.logrho < 1) {  # candidate is not symmetric
    R <- R + dnorm(logrho, can.logrho, mh.rho, log = TRUE) -
      pnorm(upper.logrho, can.logrho, mh.rho, log.p = TRUE) -
      dnorm(can.logrho, logrho, mh.rho, log = TRUE) +
      pnorm(upper.logrho, logrho, mh.rho, log.p = TRUE)
  }

  if (!fixnu) {
    R <- R + dnorm(can.lognu, lognu.m, lognu.s, log = TRUE) -
      dnorm(lognu, lognu.m, lognu.s, log = TRUE)
    if (upper.lognu < 1) {  # candidate is not symmetric
      R <- R + dnorm(lognu, can.lognu, mh.nu, log = TRUE) -
        pnorm(upper.lognu, can.lognu, mh.nu, log.p = TRUE) -
        dnorm(can.lognu, lognu, mh.nu, log = TRUE) +
        pnorm(upper.lognu, lognu, mh.nu, log.p = TRUE)
    }
  }

  if (!is.na(R)) { if (log(runif(1)) < R) {
    rho         <- can.rho
    nu          <- can.nu
    gamma       <- can.gamma
    prec        <- can.prec
    logdet.prec <- can.logdet.prec
    acc.rho     <- acc.rho + 1
    acc.nu      <- acc.nu + 1
    acc.gamma   <- acc.gamma + 1
    cur.rss     <- can.rss
  }}

  results <- list(rho = rho, nu = nu, gamma = gamma, prec = prec,
                  logdet.prec = logdet.prec, cur.rss = cur.rss,
                  att.rho = att.rho, acc.rho = acc.rho, att.nu = att.nu,
                  acc.nu = acc.nu, att.gamma = att.gamma,
                  acc.gamma = acc.gamma)
}

updateLambda1_disc_lambda <- function(x.beta, zg, y, prec, taug) {
  lll <- c(0, 0, 0)
  mmm <- c(-1, 0, 1)

  for (l in 1:length(lll)) {
    res.l  <- y - x.beta - mmm[l] * zg
    lll[l] <- -0.5 * sum(rss(prec, y = sqrt(taug) * res.l))
  }

  lambda.1 <- sample(mmm, 1, prob=exp(lll - max(lll)))
}

updateLambda2_disc_lambda <- function(lambda.a, lambda.b, z, tau) {
  nknots <- nrow(z)
  nt     <- ncol(z)
  aaa    <- lambda.a + 0.5 * nknots * nt
  bbb    <- lambda.b + 0.5 * sum(z^2 * tau)
  lambda.2 <- rgamma(1, aaa, bbb)
}

predictY_disc_lambda <- function(d11, d12, cov.model, rho, nu, gamma, res,
                                 beta, tau, taug, z, prec, lambda.1, s.pred,
                                 x.pred, knots) {
  np <- nrow(d11)
  nt <- ncol(res)
  ns <- nrow(res)
  nknots <- nrow(tau)
  yp <- matrix(NA, np, nt)

  if (cov.model == "matern") {
    s.11 <- gamma * simple.cov.sp(D = d11, sp.type = "matern",
                                  sp.par = c(1, rho),
                                  error.var = 0, smoothness = nu,
                                  finescale.var = 0)
    s.12 <- gamma * simple.cov.sp(D = d12, sp.type = "matern",
                                  sp.par = c(1, rho),
                                  error.var = 0, smoothness = nu,
                                  finescale.var = 0)
  } else {
    s.11 <- gamma * matrix(exp(-d11 / rho), np, np)
    s.12 <- gamma * matrix(exp(-d12 / rho), np, ns)
  }
  diag(s.11) <- 1
  s.12.22.inv <- s.12 %*% prec
  corp <- s.11 - s.12.22.inv %*% t(s.12)
  corp.sd.mtx <- tryCatch(chol(corp),  # only want the cholesky factor
                          error = function(e) {
                            eig.inv(corp, inv = FALSE, logdet = FALSE,
                                    mtx.sqrt = TRUE)$sd.mtx
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

predictY_cont_lambda <- function(d11, d12, cov.model, rho, nu, gamma, res,
                                 beta, tau, taug, z, prec, lambda, s.pred,
                                 x.pred, knots) {
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

# y is data matrix
# thresh.mtx is matrix of upper thresholds
# taug is site-specific precision terms
# mu is mean
# obs is which observations need to be imputed

imputeY <- function(y, taug, mu, obs, cor, gamma, thresh.mtx = NULL) {
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
        thresh.these <- Inf
        u.upper <- rep(1, length(impute.these))
      } else {
        thresh.these <- thresh.mtx[impute.these, t]
        u.upper <- pnorm(thresh.these, impute.e, impute.sd)
      }

      u.impute   <- runif(length(impute.these))
      y.impute.t <- impute.e + impute.sd * qnorm(u.impute * u.upper)
      if (sum(u.upper < 1e-6) > 0) {
        y.impute.t[(u.upper < 1e-6)] <- 0.99999 * thresh.these[(u.upper < 1e-6)]
      }
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

updateKnotsTS_disc_lambda <- function(phi, knots, g, ts, tau, z, s, min.s,
                                      max.s, x.beta, lambda.1, y, prec, att,
                                      acc, mh, update.prop = 1, att.phi,
                                      acc.phi, mh.phi) {
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- dim(knots)[1]
  avgparts <- rep(0, nt)
  taug <- zg <- matrix(NA, ns, nt)

  # recalculate knots.star at the beginning
  knots.star <- array(NA, dim = c(nknots, 2, nt))
  knots.star[, 1, ] <- transform$probit(knots[, 1, ], lower = min.s[1],
                                        upper = max.s[1])
  knots.star[, 2, ] <- transform$probit(knots[, 2, ], lower = min.s[2],
                                        upper = max.s[2])

  if (!ts) {  # will be returning these with the function results
    phi <- att.phi <- acc.phi <- mh.phi <- 0
  }

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    taug.t   <- tau[g[, t], t]
    y.t      <- y[, t]
    x.beta.t <- x.beta[, t]
    zg.t     <- z[g[, t], t]
    res.t    <- y.t - x.beta.t - lambda.1 * zg.t
    cur.lly  <- 0.5 * sum(log(taug.t)) -
      0.5 * quad.form(prec, sqrt(taug.t) * res.t)

    can.knots.star <- cur.knots.star <- knots.star[, , t]
    for (k in 1:nknots) {
      can.knots.star[k, ] <- cur.knots.star[k, ] + mh[k, t] * rnorm(2)
    }
    can.knots <- matrix(NA, nknots, 2)
    can.knots[, 1] <- transform$inv.probit(can.knots.star[, 1],
                                           lower = min.s[1],
                                           upper = max.s[1])
    can.knots[, 2] <- transform$inv.probit(can.knots.star[, 2],
                                           lower = min.s[2],
                                           upper = max.s[2])

    # recreate the partition
    can.g    <- mem(s, can.knots)
    can.taug <- tau[can.g, t]
    can.zg   <- z[can.g, t]
    can.res  <- y.t - x.beta.t - lambda.1 * can.zg
    can.lly  <- 0.5 * sum(log(can.taug)) -
      0.5 * quad.form(prec, sqrt(can.taug) * can.res)

    # remember, when not a TS, phi = 0
    if (ts & (t > 1)) {
      mean <- phi * knots.star[, , (t - 1)]
      sd   <- sqrt(1 - phi^2)
    } else {
      mean <- 0
      sd <- 1
    }

    R <- can.lly - cur.lly +
      sum(dnorm(can.knots.star, mean, sd, log = TRUE)) -
      sum(dnorm(cur.knots.star, mean, sd, log = TRUE))

    # time series also needs ot adjust R to account for next day
    if (ts & (t < nt)) {
      sd <- sqrt( 1- phi^2)
      can.mean <- phi * can.knots.star
      cur.mean <- phi * cur.knots.star
      knots.lag1 <- knots.star[, , (t + 1)]
      R <- R + sum(dnorm(knots.lag1, can.mean, sd, log = TRUE)) -
        sum(dnorm(knots.lag1, cur.mean, sd, log = TRUE))
    }

    if (!is.na(R)) { if (log(runif(1)) < R) {
      acc[, t] <- acc[, t] + 1
      knots.star[, , t] <- can.knots.star
      knots[, , t] <- can.knots
      g[, t] <- can.g
    }}

    zg[, t] <- z[g[, t], t]
    taug[, t] <- tau[g[, t], t]
  }

  if (ts) {
    phi.update <- updatePhiTS(data = knots.star, phi = phi, day.mar = 3,
                              att = att.phi, acc = acc.phi, mh = mh.phi)
    phi     <- phi.update$phi
    acc.phi <- phi.update$acc
    att.phi <- phi.update$att
  }

  results <- list(knots.star = knots.star, knots = knots, g = g, taug = taug,
                  zg = zg, acc = acc, att = att, phi = phi, acc.phi = acc.phi,
                  att.phi = att.phi)

  return(results)
}

updateKnotsTS_cont_lambda <- function(phi, knots, g, ts, tau, z, s, min.s,
                                      max.s, x.beta, lambda, y, prec, att,
                                      acc, mh, update.prop = 1, att.phi,
                                      acc.phi, mh.phi) {
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- dim(knots)[1]
  avgparts <- rep(0, nt)
  taug <- zg <- matrix(NA, ns, nt)

  # recalculate knots.star at the beginning
  knots.star <- array(NA, dim = c(nknots, 2, nt))
  knots.star[, 1, ] <- transform$probit(knots[, 1, ], lower = min.s[1],
                                        upper = max.s[1])
  knots.star[, 2, ] <- transform$probit(knots[, 2, ], lower = min.s[2],
                                        upper = max.s[2])

  if (!ts) {  # will be returning these with the function results
    phi <- att.phi <- acc.phi <- mh.phi <- 0
  }

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    taug.t   <- tau[g[, t], t]
    y.t      <- y[, t]
    x.beta.t <- x.beta[, t]
    zg.t     <- z[g[, t], t]
    res.t    <- y.t - x.beta.t - lambda * zg.t
    cur.lly  <- 0.5 * sum(log(taug.t)) -
      0.5 * quad.form(prec, sqrt(taug.t) * res.t)

    can.knots.star <- cur.knots.star <- knots.star[, , t]
    for (k in 1:nknots) {
      can.knots.star[k, ] <- cur.knots.star[k, ] + mh[k, t] * rnorm(2)
    }
    can.knots <- matrix(NA, nknots, 2)
    can.knots[, 1] <- transform$inv.probit(can.knots.star[, 1],
                                           lower = min.s[1],
                                           upper = max.s[1])
    can.knots[, 2] <- transform$inv.probit(can.knots.star[, 2],
                                           lower = min.s[2],
                                           upper = max.s[2])

    # recreate the partition
    can.g <- mem(s, can.knots)
    can.taug <- tau[can.g, t]
    can.zg <- z[can.g, t]
    can.res <- y.t - x.beta.t - lambda * can.zg
    can.lly <- 0.5 * sum(log(can.taug)) -
      0.5 * quad.form(prec, sqrt(can.taug) * can.res)

    # remember, when not a TS, phi = 0
    if (ts & (t > 1)) {
      mean <- phi * knots.star[, , (t - 1)]
      sd   <- sqrt(1 - phi^2)
    } else {
      mean <- 0
      sd <- 1
    }

    R <- can.lly - cur.lly +
      sum(dnorm(can.knots.star, mean, sd, log = TRUE)) -
      sum(dnorm(cur.knots.star, mean, sd, log = TRUE))

    # time series also needs to adjust R to account for next day
    if (ts & (t < nt)) {
      sd <- sqrt(1 - phi^2)
      can.mean <- phi * can.knots.star
      cur.mean <- phi * cur.knots.star
      knots.lag1 <- knots.star[, , (t + 1)]
      R <- R + sum(dnorm(knots.lag1, can.mean, sd, log = TRUE)) -
        sum(dnorm(knots.lag1, cur.mean, sd, log = TRUE))
    }

    if (!is.na(R)) { if (log(runif(1)) < R) {
      acc[, t] <- acc[, t] + 1
      knots.star[, , t] <- can.knots.star
      knots[, , t] <- can.knots
      g[, t] <- can.g
    }}

    zg[, t] <- z[g[, t], t]
    taug[, t] <- tau[g[, t], t]
  }

  if (ts) {
    phi.update <- updatePhiTS(data = knots.star, phi = phi, day.mar = 3,
                              att = att.phi, acc = acc.phi, mh = mh.phi)
    phi     <- phi.update$phi
    acc.phi <- phi.update$acc
    att.phi <- phi.update$att
  }

  results <- list(knots.star = knots.star, knots = knots, g = g, taug = taug,
                  zg = zg, acc = acc, att = att, phi = phi, acc.phi = acc.phi,
                  att.phi = att.phi)

  return(results)
}

# data is the data for the time series
# phi is the parameter for the time series
#   data_t ~ N(phi * data_{t-1}, sqrt(1 - phi^2))
# att.phi is the current number of attempts
# acc.phi is the current number of acceptances
# mh.phi is the candidate standard deviation
# day.mar is the margin in the data that represents time

# MAKE SURE THAT TIME IS LAST MARGIN OF DATA
updatePhiTS <- function(data, phi, day.mar, att, acc, mh) {
  att <- att + 1
  nt <- dim(data)[day.mar]
  if (day.mar == 2) {
    data.up1  <- data[, -1]
    data.lag1 <- data[, -nt]
  } else if (day.mar == 3) {
    data.up1  <- data[, , -1]
    data.lag1 <- data[, , -nt]
  }

  cur.mean <- phi * data.lag1
  cur.sd   <- sqrt(1 - phi^2)

  cur.phi.star <- transform$probit(x = phi, lower = -1, upper = 1)
  can.phi.star <- rnorm(1, cur.phi.star, mh)
  can.phi      <- transform$inv.probit(x = can.phi.star, lower = -1, upper = 1)
  can.mean     <- can.phi * data.lag1
  can.sd       <- sqrt(1 - can.phi^2)

  # likelihood of data impacted by phi does not include day 1
  R <- sum(dnorm(data.up1, can.mean, can.sd, log = TRUE)) -
    sum(dnorm(data.up1, cur.mean, cur.sd, log = TRUE)) +
    dnorm(can.phi.star, log = TRUE) - dnorm(cur.phi.star, log = TRUE)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc <- acc + 1
    phi <- can.phi
  }}

  results <- list(phi = phi, att = att, acc = acc)
}