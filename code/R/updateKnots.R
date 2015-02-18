updateKnotsTS <- function(phi, knots, g, ts, tau, z, s, min.s, max.s, x.beta,
                          lambda.1, y, prec, att, acc, mh, update.prop=1,
                          att.phi, acc.phi, mh.phi) {
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- dim(knots)[1]
  avgparts <- rep(0, nt)
  taug <- zg <- matrix(NA, ns, nt)

  # recalculate knots.star at the beginning
  knots.star <- array(NA, dim=c(nknots, 2, nt))
  knots.star[, 1, ] <- transform$probit(knots[, 1, ], lower=min.s[1],
                                        upper=max.s[1])
  knots.star[, 2, ] <- transform$probit(knots[, 2, ], lower=min.s[2],
                                        upper=max.s[2])

  if (!ts) {  # will be returning these with the function results
    phi <- att.phi <- acc.phi <- mh.phi <- 0
  }

  for (t in 1:nt) {
    taug.t   <- tau[g[, t], t]
    y.t      <- y[, t]
    x.beta.t <- x.beta[, t]
    zg.t     <- z[g[, t], t]
    res.t    <- y.t - x.beta.t - lambda.1 * zg.t
    cur.lly  <- 0.5 * sum(log(taug.t)) -
                0.5 * quad.form(prec, sqrt(taug.t) * res.t)

    att[, t] <- att[, t] + 1
    can.knots.star <- cur.knots.star <- knots.star[, , t]
    for (k in 1:nknots) {
      can.knots.star[k, ] <- cur.knots.star[k, ] + mh[k, t] * rnorm(2)
    }
    can.knots <- matrix(NA, nknots, 2)
    can.knots[, 1] <- transform$inv.probit(can.knots.star[, 1], lower=min.s[1],
                                           upper=max.s[1])
    can.knots[, 2] <- transform$inv.probit(can.knots.star[, 2], lower=min.s[2],
                                           upper=max.s[2])

    # recreate the partition
    can.g <- mem(s, can.knots)
    can.taug <- tau[can.g, t]
    can.zg <- z[can.g, t]
    can.res <- y.t - x.beta.t - lambda.1 * can.zg
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
         sum(dnorm(can.knots.star, mean, sd, log=T)) -
         sum(dnorm(cur.knots.star, mean, sd, log=T))

    # time series also needs ot adjust R to account for next day
    if (ts & (t < nt)) {
      sd <- sqrt( 1- phi^2)
      can.mean <- phi * can.knots.star
      cur.mean <- phi * cur.knots.star
      knots.lag1 <- knots.star[, , (t + 1)]
      R <- R + sum(dnorm(knots.lag1, can.mean, sd, log=T)) -
               sum(dnorm(knots.lag1, cur.mean, sd, log=T))
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

    # for (k in 1:nknots) {
    #   sample <- rbinom(1, 1, update.prop)
    #   if (sample) {
    #     att[k, t] <- att[k, t] + 1

    #     cur.knots.star      <- knots.star[, , t]
    #     can.knots.star      <- cur.knots.star
    #     can.knots.star[k, ] <- cur.knots.star[k, ] + mh[k, t] * rnorm(2)
    #     can.knots           <- matrix(NA, nknots, 2)
    #     can.knots[, 1]      <- transform$inv.probit(can.knots.star[, 1],
    #                                                 lower=min.s[1],
    #                                                 upper=max.s[1])
    #     can.knots[, 2]      <- transform$inv.probit(can.knots.star[, 2],
    #                                                 lower=min.s[2],
    #                                                 upper=max.s[2])

    #     # recreate the partition
    #     can.g    <- mem(s, can.knots)
    #     can.taug <- tau[can.g, t]
    #     can.zg   <- z[can.g, t]
    #     can.res  <- y.t - x.beta.t - lambda.1 * can.zg
    #     can.lly  <- 0.5 * sum(log(can.taug)) -
    #                 0.5 * quad.form(prec, sqrt(can.taug) * can.res)

    #     # remember, when not a TS, phi = 0
    #     if (ts & (t > 1)) {
    #       mean <- phi * knots.star[k, , (t - 1)]
    #       sd   <- sqrt(1 - phi^2)
    #     } else {
    #       mean <- 0
    #       sd   <- 1
    #     }

    #     R <- can.lly - cur.lly +
    #          sum(dnorm(can.knots.star[k, ], mean, sd, log=T)) -
    #          sum(dnorm(cur.knots.star[k, ], mean, sd, log=T))

    #     # time series also needs to adjust R to account for next day
    #     if (ts & (t < nt)) {
    #       sd <- sqrt(1 - phi^2)
    #       can.mean   <- phi * can.knots.star[k, ]
    #       cur.mean   <- phi * cur.knots.star[k, ]
    #       knots.lag1 <- knots.star[k, , (t + 1)]
    #       R <- R + sum(dnorm(knots.lag1, can.mean, sd, log=T)) -
    #                sum(dnorm(knots.lag1, cur.mean, sd, log=T))
    #     }

    #     if (!is.na(R)) { if (log(runif(1)) < R) {
    #       acc[k, t]          <- acc[k, t] + 1
    #       knots.star[k, , t] <- can.knots.star[k, ]
    #       knots[k, , t]      <- can.knots[k, ]
    #       g[, t]             <- can.g
    #       taug.t             <- can.taug
    #       zg.t               <- can.zg
    #       cur.lly            <- can.lly
    #     }}
    #   }
    #   zg[, t]   <- zg.t
    #   taug[, t] <- taug.t
    # }

  # }

  if (ts) {
    phi.update <- updatePhiTS(data=knots.star, phi=phi, day.mar=3,
                              att=att.phi, acc=acc.phi, mh=mh.phi)
    phi     <- phi.update$phi
    acc.phi <- phi.update$acc
    att.phi <- phi.update$att
  }

  results <- list(knots.star=knots.star, knots=knots, g=g, taug=taug, zg=zg,
                  acc=acc, att=att, phi=phi, acc.phi=acc.phi, att.phi=att.phi)

  return(results)
}

updateKnotsTS1 <- function(phi, knots, g, ts, tau, z, s, min.s, max.s, x.beta,
                           lambda, y, prec, att, acc, mh, update.prop=1,
                           att.phi, acc.phi, mh.phi) {
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- dim(knots)[1]
  avgparts <- rep(0, nt)
  taug <- zg <- matrix(NA, ns, nt)

  # recalculate knots.star at the beginning
  knots.star <- array(NA, dim=c(nknots, 2, nt))
  knots.star[, 1, ] <- transform$probit(knots[, 1, ], lower=min.s[1],
                                        upper=max.s[1])
  knots.star[, 2, ] <- transform$probit(knots[, 2, ], lower=min.s[2],
                                        upper=max.s[2])

  if (!ts) {  # will be returning these with the function results
    phi <- att.phi <- acc.phi <- mh.phi <- 0
  }

  for (t in 1:nt) {
    taug.t   <- tau[g[, t], t]
    y.t      <- y[, t]
    x.beta.t <- x.beta[, t]
    zg.t     <- z[g[, t], t]
    res.t    <- y.t - x.beta.t - lambda * zg.t
    cur.lly  <- 0.5 * sum(log(taug.t)) -
                0.5 * quad.form(prec, sqrt(taug.t) * res.t)

    att[, t] <- att[, t] + 1
    can.knots.star <- cur.knots.star <- knots.star[, , t]
    for (k in 1:nknots) {
      can.knots.star[k, ] <- cur.knots.star[k, ] + mh[k, t] * rnorm(2)
    }
    can.knots <- matrix(NA, nknots, 2)
    can.knots[, 1] <- transform$inv.probit(can.knots.star[, 1], lower=min.s[1],
                                           upper=max.s[1])
    can.knots[, 2] <- transform$inv.probit(can.knots.star[, 2], lower=min.s[2],
                                           upper=max.s[2])

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
         sum(dnorm(can.knots.star, mean, sd, log=T)) -
         sum(dnorm(cur.knots.star, mean, sd, log=T))

    # time series also needs to adjust R to account for next day
    if (ts & (t < nt)) {
      sd <- sqrt( 1- phi^2)
      can.mean <- phi * can.knots.star
      cur.mean <- phi * cur.knots.star
      knots.lag1 <- knots.star[, , (t + 1)]
      R <- R + sum(dnorm(knots.lag1, can.mean, sd, log=T)) -
               sum(dnorm(knots.lag1, cur.mean, sd, log=T))
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
    phi.update <- updatePhiTS(data=knots.star, phi=phi, day.mar=3,
                              att=att.phi, acc=acc.phi, mh=mh.phi)
    phi     <- phi.update$phi
    acc.phi <- phi.update$acc
    att.phi <- phi.update$att
  }

  results <- list(knots.star=knots.star, knots=knots, g=g, taug=taug, zg=zg,
                  acc=acc, att=att, phi=phi, acc.phi=acc.phi, att.phi=att.phi)

  return(results)
}


# att[, t] <- att[, t] + 1
#     can.knots.star <- knots.star[, , t]
#     for (k in 1:nknots) {
#       can.knots.star[k, ] <- can.knots.star[k, ] + mh[k, t] * rnorm(2)
#     }
#     can.knots <- matrix(NA, nknots, 2)
#     can.knots[, 1] <- transform$inv.probit(can.knots.star[, 1], lower=min.s[1],
#                                            upper=max.s[1])
#     can.knots[, 2] <- transform$inv.probit(can.knots.star[, 2], lower=min.s[2],
#                                            upper=max.s[2])

#     # recreate the partition
#     can.g    <- mem(s, can.knots)
#     can.taug <- tau[can.g, t]
#     can.zg   <- z[can.g, t]
#     can.res  <- y.t - x.beta.t - lambda.1 * can.zg
#     can.lly  <- 0.5 * sum(log(can.taug)) -
#                 0.5 * quad.form(prec, sqrt(can.taug) * can.res)

#     # remember, when not a TS, phi = 0
#     if (ts & (t > 1)) {
#       mean <- phi * knots.star[, , (t - 1)]
#       sd   <- sqrt(1 - phi^2)
#     } else {
#       mean <- 0
#       sd   <- 1
#     }

#     R <- can.lly - cur.lly +
#          sum(dnorm(can.knots.star, mean, sd, log=T)) -
#          sum(dnorm(knots.star[, , t], mean, sd, log=T))

#     # time series also needs to adjust R to account for next day
#     if (ts & (t < nt)) {
#       sd       <- sqrt(1 - phi^2)
#       can.mean <- phi * can.knots.star
#       cur.mean <- phi * knots.star[, , t]
#       knots.lag1 <- knots.star[, , (t + 1)]
#       R <- R + sum(dnorm(knots.lag1, can.mean, sd, log=T)) -
#                sum(dnorm(knots.lag1, cur.mean, sd, log=T))
#     }

#     if (!is.na(R)) { if (log(runif(1)) < R) {
#       acc[, t]          <- acc[, t] + 1
#       knots.star[, , t] <- can.knots.star
#       knots[, , t]      <- can.knots
#       g[, t]            <- can.g
#       taug.t            <- can.taug
#       zg.t              <- can.zg
#       # cur.lly            <- can.lly
#     }}

#     zg[, t]   <- zg.t
#     taug[, t] <- taug.t