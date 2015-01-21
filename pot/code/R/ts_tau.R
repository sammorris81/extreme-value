updateTauTS <- function(phi, tau, taug, g, res, nparts.tau, prec, z,
                        tau.alpha, tau.beta, skew,
                        att, acc, mh, att.tau, acc.tau,
                        att.phi, acc.phi, mh.phi) {
  ns <- nrow(res)
  nt <- ncol(res)
  nknots <- nrow(tau)

  tau.star <- qnorm(pgamma(tau, tau.alpha, tau.beta))

  if (nknots == 1) {  # seems to come back ok
    for (t in 1:nt) {
      att[1, t] <- att[1, t] + 1
      res.t     <- res[, t]
      cur.lly   <- 0.5 * ns * log(tau[1, t]) -
                   0.5 * tau[1, t] * quad.form(prec, res.t)

      can.tau.star  <- rnorm(1, tau.star[1, t], mh[1, t])

      # transform back to R+
      can.tau <- qgamma(pnorm(can.tau.star), tau.alpha, tau.beta)

      can.lly <- 0.5 * ns * log(can.tau) -
                 0.5 * can.tau * quad.form(prec, res.t)

      if (skew) {
        cur.llz <- 0.5 * log(tau[1, t]) - 0.5 * tau[1, t] * z[1, t]^2
        can.llz <- 0.5 * log(can.tau) - 0.5 * can.tau * z[1, t]^2
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
      R <- R + dnorm(can.tau.star, mean, sd, log=T) -
               dnorm(tau.star[1, t], mean, sd, log=T)

      # adjust R to account for next day
      if (t < nt) {
        sd            <- sqrt(1 - phi^2)
        can.mean      <- phi * can.tau.star
        cur.mean      <- phi * tau.star[1, t]
        tau.star.lag1 <- tau.star[1, t + 1]

        R <- R + dnorm(tau.star.lag1, can.mean, sd, log=T) -
                 dnorm(tau.star.lag1, cur.mean, sd, log=T)
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
      res.t   <- res[, t]
      cur.lly <- 0.5 * sum(log(taug[, t])) -
                 0.5 * quad.form(prec, sqrt(taug[, t]) * res.t)

      for (k in 1:nknots) {
        these  <- which(g[, t] == k)
        nparts <- length(these)
        nparts.tau[k, t] <- nparts

        can.tau.star <- tau.star[, t]
        can.tau.star[k] <- rnorm(1, tau.star[k, t], mh[k, t])

        # transform back to R+
        can.tau  <- qgamma(pnorm(can.tau.star), tau.alpha, tau.beta)
        can.taug <- can.tau[g[, t]]

        can.lly <- 0.5 * sum(log(can.taug)) -
                   0.5 * quad.form(prec, sqrt(can.taug) * res.t)

        if (skew) {
          cur.llz <- 0.5 * log(tau[k, t]) - 0.5 * tau[k, t] * z[k, t]^2
          can.llz <- 0.5 * log(can.tau[k]) - 0.5 * can.tau[k] * z[k, t]^2
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

        R <- R + dnorm(can.tau.star, mean, sd, log=T) -
                 dnorm(tau.star[k, t], mean, sd, log=T)

        if (t < nt) {
          sd            <- sqrt(1 - phi^2)
          can.mean      <- phi * can.tau.star
          cur.mean      <- phi * tau.star[k, t]
          tau.star.lag1 <- tau.star[k, t + 1]
          R <- R + dnorm(tau.star.lag1, can.mean, sd, log=T) =
                   dnorm(tau.star.lag1, cur.mean, sd, log=T)
        }

        if (!is.na(R)) { if (log(runif(1)) < R) {
          acc[k, t]      <- acc[k, t] + 1
          acc.tau[k, t]  <- acc.tau[k, t] + 1
          tau.star[k, t] <- can.tau.star[k]
          tau[k, t]      <- can.tau[k]
          taug[these, t] <- can.tau[k]
          cur.lly        <- can.lly
        }}
      }  # end k
    }  # end t
  }  # fi nknots > 1

  att.phi       <- att.phi + 1
  tau.star.lag1 <- tau.star[, -nt]             # don't need the last day
  cur.mean      <- phi * tau.star.lag1
  cur.sd        <- sqrt(1 - phi^2)

  phi.con       <- qnorm((phi + 1) / 2)        # transform to R
  can.phi.con   <- rnorm(1, phi.con, mh.phi)   # draw candidate
  can.phi       <- 2 * pnorm(can.phi.con) - 1  # transform back to (-1, 1)
  can.mean      <- can.phi * tau.star.lag1
  can.sd        <- sqrt(1 - can.phi^2)

  # likelihood impacted by phi.tau does not include first day
  R <- sum(dnorm(tau.star[, -1], can.mean, can.sd, log=T)) -
       sum(dnorm(tau.star[, -1], cur.mean, cur.sd, log=T)) +
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi     <- can.phi
  } }

  results <- list(tau=tau, taug=taug, phi=phi, acc=acc, att=att,
                  acc.tau=acc.tau, att.tau=att.tau,
                  acc.phi=acc.phi, att.phi=att.phi)

  # results <- list(phi=phi, acc.phi=acc.phi, att.phi=att.phi)
}