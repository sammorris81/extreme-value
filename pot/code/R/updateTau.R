updateTauGaus <- function(res, prec, tau.alpha, tau.beta) {
  ns <- nrow(res)
  nt <- ncol(res)

  this.rss <- sum(rss(prec=prec, y=res))
  aaa <- tau.alpha + 0.5 * ns * nt
  bbb <- tau.beta + 0.5 * this.rss
  tau <- rgamma(1, aaa, bbb)

  return(tau)
}

updateTau <- function(tau, taug, g, res, nparts.tau, prec, z, lambda.2,
                      tau.alpha, tau.beta, skew,
                      att, acc, mh) {
  ns <- nrow(res)
  nt <- ncol(res)
  nknots <- nrow(tau)

  if (nknots == 1) {
    for (t in 1:nt) {
      res.t <- res[, t]
      rss.t <- quad.form(prec, res.t)

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
      res.t <- res[, t]
      cur.lly <- 0.5 * sum(log(taug[, t])) -
                 0.5 * quad.form(prec, sqrt(taug[, t]) * res.t)

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
          bbb <- tau.beta + 0.5 * quad.form(prec, res.t)

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
          bbb <- 0.5 * quad.form(prec[these, these], res.t[these])

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
                     0.5 * quad.form(prec, sqrt(can.taug) * res.t)

          if (skew) {
            cur.llz <- 0.5 * log(tau[k, t]) -
                       0.5 * tau[k, t] * z[k, t]^2 * lambda.2
            can.llz <- 0.5 * log(can.tau[k]) -
                       0.5 * can.tau[k] * z[k, t]^2 * lambda.2
          } else {
            cur.llz <- can.llz <- 0
          }

          R <- can.lly - cur.lly + can.llz - cur.llz +
               dgamma(can.tau[k], tau.alpha, tau.beta, log=T) -
               dgamma(tau[k, t], tau.alpha, tau.beta, log=T) +
               dgamma(tau[k, t], aaa, bbb, log=T) -
               dgamma(can.tau[k], aaa, bbb, log=T)


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

  results <- list(tau=tau, taug=taug, acc=acc, att=att)

}

updateTauTS <- function(phi, tau, taug, g, res, nparts.tau, prec, z, lambda.2,
                        tau.alpha, tau.beta, skew, att, acc, mh,
                        att.phi, acc.phi, mh.phi) {
  ns <- nrow(res)
  nt <- ncol(res)
  nknots <- nrow(tau)

  tau.star <- gamma.cop(tau, tau.alpha, tau.beta)

  if (nknots == 1) {  # seems to come back ok
    for (t in 1:nt) {
      att[1, t] <- att[1, t] + 1
      res.t     <- res[, t]
      cur.lly   <- 0.5 * ns * log(tau[1, t]) -
                   0.5 * tau[1, t] * quad.form(prec, res.t)

      can.tau.star  <- rnorm(1, tau.star[1, t], mh[1, t])

      # transform back to R+
      can.tau <- gamma.invcop(can.tau.star, tau.alpha, tau.beta)

      can.lly <- 0.5 * ns * log(can.tau) -
                 0.5 * can.tau * quad.form(prec, res.t)

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
        att[k, t] <- att[k, t] + 1
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

        R <- R + dnorm(can.tau.star[k], mean, sd, log=T) -
                 dnorm(tau.star[k, t], mean, sd, log=T)

        if (t < nt) {
          sd            <- sqrt(1 - phi^2)
          can.mean      <- phi * can.tau.star[k]
          cur.mean      <- phi * tau.star[k, t]
          tau.star.lag1 <- tau.star[k, t + 1]
          R <- R + dnorm(tau.star.lag1, can.mean, sd, log=T) -
                   dnorm(tau.star.lag1, cur.mean, sd, log=T)
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

  phi.update <- updatePhiTS(data=tau.star, phi=phi, day.mar=2,
                            att=att.phi, acc=acc.phi, mh=mh.phi)
  phi     <- phi.update$phi
  acc.phi <- phi.update$acc
  att.phi <- phi.update$att

  results <- list(tau=tau, taug=taug, phi=phi, acc=acc, att=att,
                  acc.phi=acc.phi, att.phi=att.phi)

  # results <- list(phi=phi, acc.phi=acc.phi, att.phi=att.phi)
}

updateTauAlpha <- function(tau, tau.beta) {
  # using unique so we can use the same function for all possible setups
  tau <- unique(tau)
  lll <- mmm <- seq(0.1, 10, 0.1)
  for (l in 1:length(lll)) {
    lll[l] <- sum(dgamma(tau, mmm[l], tau.beta, log=T))
  }

  tau.alpha <- sample(mmm, 1, prob=exp(lll - max(lll)))
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

updateTauOld <- function(y, mu, tau, taug, g, nparts.tau, prec, z, lambda.2,
                      tau.alpha, tau.beta, skew,
                      att, acc, mh, att.tau, acc.tau, mh.tau) {
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- nrow(tau)
  res <- y - mu

  if (nknots == 1) {
    for (t in 1:nt) {
      res.t <- res[, t]
      rss.t <- quad.form(prec, res.t)

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
      res.t <- res[, t]
      cur.lly <- 0.5 * sum(log(taug[, t])) -
                 0.5 * quad.form(prec, sqrt(taug[, t]) * res.t)

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
          bbb <- tau.beta + 0.5 * quad.form(prec, res.t)

          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2 * lambda.2
          }

          tau[k, t] <- rgamma(1, aaa, bbb)
          taug[, t] <- tau[k, t]

        } else {  # nparts > 0
          att[nparts + 1] <- att[nparts + 1] + 1
          att.tau[k, t] <- att.tau[k, t] + 1

          can.logtau    <- log(tau[, t])
          can.logtau[k] <- rnorm(1, log(tau[k, t]), mh.tau[k, t])
          can.tau  <- exp(can.logtau)
          can.taug <- taug[, t]
          can.taug[these] <- can.tau[k]

          can.lly <- 0.5 * sum(log(can.taug)) -
                     0.5 * quad.form(prec, sqrt(can.taug) * res.t)

          if (skew) {
            cur.llz <- 0.5 * log(tau[k, t]) -
                       0.5 * tau[k, t] * z[k, t]^2 * lambda.2
            can.llz <- 0.5 * log(can.tau[k]) -
                       0.5 * can.tau[k] * z[k, t]^2 * lambda.2
          } else {
            cur.llz <- can.llz <- 0
          }

          R <- can.lly - cur.lly + can.llz - cur.llz +
               dgamma(can.tau[k], tau.alpha, tau.beta, log=T) -
               dgamma(tau[k, t], tau.alpha, tau.beta, log=T)

          if (!is.na(R)) { if (log(runif(1)) < R) {
            acc[nparts + 1] <- acc[nparts + 1] + 1
            acc.tau[k, t]   <- acc.tau[k, t] + 1
            tau[k, t]       <- can.tau[k]
            taug[these, t]  <- can.tau[k]
            cur.lly         <- can.lly
          }}
        }  # fi nparts
      }  # end k
    }  # end t
  }  # fi nknots > 1

  results <- list(tau=tau, taug=taug, acc=acc, att=att,
                  acc.tau=acc.tau, att.tau=att.tau)

}