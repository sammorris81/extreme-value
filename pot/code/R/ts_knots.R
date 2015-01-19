updateKnotsTS <- function(phi, knots, knots.con, g, ts, tau, taug,
                          z, zg, s, s.range, x.beta, lambda, y, res, prec,
                          att, acc, mh, att.phi, acc.phi, mh.phi) {
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- dim(knots)[1]
  avgparts <- rep(0, nt)

  if (!ts) {  # will be returning these with the function results
    phi <- att.phi <- acc.phi <- mh.phi <- 0
  }

  for (t in 1:nt) {
    att[1] <- att[1] + 1

    cur.knots.con  <- knots.con[, , t]
    can.knots.con  <- cur.knots.con + mh[1] * rnorm(2 * nknots)
    can.knots      <- pnorm(can.knots.con)
    can.knots[, 1] <- can.knots[, 1] * s.range + min(s[, 1])
    can.knots[, 2] <- can.knots[, 2] * s.range + min(s[, 2])

    # recreate the partition
    can.g          <- mem(s, can.knots)
    can.taug       <- tau[can.g, t]
    can.zg         <- z[can.g, t]
    can.res        <- y[, t] - x.beta[, t] - lambda * can.zg

    R <- -0.5 * quad.form(prec, sqrt(can.taug) * can.res) +
          0.5 * quad.form(prec, sqrt(taug[, t]) * res[, t]) +
          0.5 * sum(log(can.taug)) - 0.5 * sum(log(taug[, t]))

    # remember, when not a TS, phi = 0
    if (t > 1) {
      mean <- phi * knots.con[, , (t - 1)]
      sd   <- sqrt(1 - phi^2)
    } else {
      mean <- 0
      sd   <- 1
    }

    R <- R + sum(dnorm(can.knots.con, mean, sd, log=T)) -
             sum(dnorm(cur.knots.con, mean, sd, log=T))

    # time series also needs to adjust R to account for next day
    if (ts & (t < nt)) {
      sd <- sqrt(1 - phi^2)
      can.mean   <- phi * can.knots.con
      cur.mean   <- phi * cur.knots.con
      knots.lag1 <- knots.con[, , (t + 1)]
      R <- R + sum(dnorm(knots.lag1, can.mean, sd, log=T)) -
               sum(dnorm(knots.lag1, cur.mean, sd, log=T))
    }

    if (!is.na(R)) { if (log(runif(1)) < R) {
      knots.con[, , t] <- can.knots.con
      knots[, , t]     <- can.knots
      g[, t]           <- can.g
      taug[, t]        <- can.taug
      zg[, t]          <- can.zg
      acc[1]           <- acc[1] + 1
    }}
  }

  if (ts) {
    att.phi        <- att.phi + 1
    knots.con.lag1 <- knots.con[, , -nt]  # don't need the last day
    cur.mean       <- phi * knots.con.lag1
    cur.sd         <- sqrt(1 - phi^2)

    phi.con     <- qnorm((phi + 1) / 2)        # transform to R
    can.phi.con <- rnorm(1, phi.con, mh.phi)   # draw candidate
    can.phi     <- 2 * pnorm(can.phi.con) - 1  # transform back to (-1, 1)
    can.mean    <- can.phi * knots.con.lag1
    can.sd      <- sqrt(1 - can.phi^2)

    # likelihood impacted by phi.w does not include first day
    R <- sum(dnorm(knots.con[, , -1], can.mean, can.sd, log=T)) -
         sum(dnorm(knots.con[, , -1], cur.mean, cur.sd, log=T)) +
         dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)

    if (!is.na(R)) { if (log(runif(1)) < R) {
      acc.phi <- acc.phi + 1
      phi     <- can.phi
    }}
  }

  results <- list(knots.con=knots.con, knots=knots, g=g, taug=taug, zg=zg,
                  acc=acc, att=att, phi=phi, acc.phi=acc.phi, att.phi=att.phi)

  return(results)
}