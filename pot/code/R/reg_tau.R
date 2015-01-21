updateTau <- function(tau, taug, g, res, nparts.tau, prec, z,
                      tau.alpha, tau.beta, skew,
                      att, acc, mh, att.tau, acc.tau) {
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
        bbb <- bbb + 0.5 * z[1, t]^2
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

        if (nparts == 0) {  # tau update is conjugate
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
        } else {  # nparts > 0
          att[nparts + 1] <- att[nparts + 1] + 1
          att.tau[k, t] <- att.tau[k, t] + 1

          aaa <- tau.alpha + 0.5 * nparts
          bbb <- tau.beta + 0.5 * quad.form(prec[these, these], res.t[these])

          if (skew) {
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[k, t]^2
          }

          # the posterior is conjugate when nparts = ns
          # as nparts -> 1, want a wider candidate
          aaa <- aaa / mh[nparts + 1]
          bbb <- bbb / mh[nparts + 1]
          if (is.nan(bbb)) {
            print(tau)
            print(z)
            print(aaa)
            print(y[1:5, ])
            print(res[1:5, ])
            print(mu[1:5, ])
            print(beta)
            print(quad.form(prec[these, these], res.t[these]))
          }
          can.tau    <- tau[, t]
          can.tau[k] <- rgamma(1, aaa, bbb)
          if (can.tau[k] < 1e-6) {
            can.tau[k] <- 1e-6
          }
          can.taug   <- can.tau[g[, t]]

          can.lly <- 0.5 * sum(log(can.taug)) -
                     0.5 * quad.form(prec, sqrt(can.taug) * res.t)

          if (skew) {
            cur.llz <- 0.5 * log(tau[k, t]) - 0.5 * tau[k, t] * z[k, t]^2
            can.llz <- 0.5 * log(can.tau[k]) - 0.5 * can.tau[k] * z[k, t]^2
          } else {
            cur.llz <- can.llz <- 0
          }

          R <- can.lly - cur.lly + can.llz - cur.llz +
               tryCatch({ dgamma(tau[k, t], aaa, bbb, log=T) },
                  warning = function(w) {
                    print(paste("knot", k, ", day", t))
                    print(paste("aaa =", aaa))
                    print(paste("bbb =", bbb))
                    print(paste("tau[k, t] =", tau[k, t]))
                    print(paste("g[, t] =", g[, t]))
                    print(paste("mh.tau[nparts + 1] =", mh.tau[nparts + 1]))
                    print(paste("tau.beta =", tau.beta))
                  }) -
               dgamma(can.tau[k], aaa, bbb, log=T)

          R <- R + dgamma(can.tau[k], tau.alpha, tau.beta, log=T) -
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