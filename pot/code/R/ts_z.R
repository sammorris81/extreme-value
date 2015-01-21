updateZTS <- function(z, zg, y, lambda, x.beta,
                      phi, tau, taug, g, prec,
                      acc, att, mh, acc.phi, att.phi, mh.phi) {
  nt <- ncol(y)
  nknots <- nrow(z)

  # transform via copula to normal
  sd     <- 1 / sqrt(tau)
  z.star <- qnorm(2 * pnorm(z, 0, sd) - 1)

  for (t in 1:nt) {
  	cur.res   <- y[, t] - x.beta[, t] + lambda * zg[, t]
    cur.rss   <- quad.form(prec, sqrt(taug[, t]) * cur.res)

    for (k in 1:nknots) {
      att[k, t]     <- att[k, t] + 1
      can.z.star    <- z.star[, t] # will be nknots long
      can.z.star[k] <- rnorm(1, z.star[k, t], mh[k, t])
      can.z         <- qnorm((pnorm(can.z.star[k]) + 1) / 2, 0, sd[k, t])
      can.zg        <- can.z[g[, t]]  # will be ns long

      can.res <- y[, t] - x.beta[, t] + lambda * can.zg
      can.rss <- quad.form(prec, sqrt(taug[, t]) * can.res)

      # prior
      if (t > 1) {
        mean <- phi * z.star[k, (t - 1)]
        sd   <- sqrt(1 - phi^2)
      } else {
        mean <- 0
        sd   <- 1
      }

      R <- -0.5 * sum(can.rss - cur.rss) +
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
        z.star[k, t] <- can.z.star[k]
        zg[, t]      <- can.zg
        cur.res      <- can.res
        cur.rss      <- can.rss
      } }
    }
  }

  # phi.z
  att.phi     <- att.phi + 1
  z.lag1      <- z.star[, -nt]               # don't need the last day
  cur.mean    <- phi * z.lag1
  cur.sd      <- sqrt(1 - phi^2)

  phi.con     <- qnorm((phi + 1) / 2)        # transform to R
  can.phi.con <- rnorm(1, phi.con, mh.phi)   # draw candidate
  can.phi     <- 2 * pnorm(can.phi.con) - 1  # transform back to (-1, 1)
  can.mean    <- can.phi * z.lag1
  can.sd      <- sqrt(1 - can.phi^2)

  # likelihood impacted by phi.z does not include first day
  R <- sum(dnorm(z.star[, -1], can.mean, can.sd, log=T)) -
       sum(dnorm(z.star[, -1], cur.mean, cur.sd, log=T)) +
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi     <- can.phi
  } }

  results <- list(z.star=z.star, acc=acc, att=att, zg=zg,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi)

  # results <- list(phi=phi, att.phi=att.phi, acc.phi=acc.phi)

  return(results)
}