updateZTS <- function(z, zg, y, lambda, x.beta,
                      phi, tau, taug, g, prec,
                      acc, att, mh, acc.phi, att.phi, mh.phi) {
  nt <- ncol(z)
  nknots <- nrow(z)

  # transform via copula to normal
  z.star <- qnorm(2 * pnorm(z * sqrt(tau)) - 1)  # nknots x nt

  for (t in 1:nt) {
  	taug.t  <- sqrt(taug[, t])
    mu.t    <- x.beta[, t] + lambda * zg[, t]
    cur.res <- y[, t] - mu.t
    cur.rss <- quad.form(prec, taug.t * cur.res)

    for (k in 1:nknots) {
      att[k, t]     <- att[k, t] + 1
      these         <- which(g[, t] == k)
      can.z.star    <- rnorm(1, z.star[k, t], mh[k, t])
      can.z         <- z[, t]
      can.z[k]      <- qnorm(0.5 * pnorm(can.z.star) + 0.5) / sqrt(tau[k, t])
      can.zg        <- can.z[g[, t]]  # ns long

      can.mu.t <- x.beta[, t] + lambda * can.zg
      can.res <- y[, t] - can.mu.t
      can.rss <- quad.form(prec, taug.t * can.res)

      # prior
      if (t > 1) {
        mean <- phi * z.star[k, (t - 1)]
        sd   <- sqrt(1 - phi^2)
      } else {
        mean <- 0
        sd   <- 1
      }

      R <- -0.5 * sum(can.rss - cur.rss) +
            dnorm(can.z.star, mean, sd, log=T) -
            dnorm(z.star[k, t], mean, sd, log=T)

      if (t < nt) {
      	sd <- sqrt(1 - phi^2)
        can.mean <- phi * can.z.star
        cur.mean <- phi * z.star[k, t]
        z.star.lag1 <- z.star[k, t + 1]
      	R <- R + dnorm(z.star.lag1, can.mean, sd, log=T) -
                 dnorm(z.star.lag1, cur.mean, sd, log=T)
      }

      if (!is.na(R)) { if (log(runif(1)) < R) {
        acc[k, t]    <- acc[k, t] + 1
        z[k, t]      <- can.z[k]
        zg[these, t] <- can.zg[these]
        z.star[k, t] <- can.z.star
        cur.res      <- can.res
        cur.rss      <- can.rss
      } }
    }
  }

  z.star <- qnorm(2 * pnorm(z * sqrt(tau)) - 1)
  # phi.z
  att.phi     <- att.phi + 1
  z.star.lag1 <- z.star[, -nt]               # don't need the last day
  cur.mean    <- phi * z.star.lag1
  cur.sd      <- sqrt(1 - phi^2)

  phi.con     <- qnorm((phi + 1) / 2)        # transform to R
  can.phi.con <- rnorm(1, phi.con, mh.phi)   # draw candidate
  can.phi     <- 2 * pnorm(can.phi.con) - 1  # transform back to (-1, 1)
  can.mean    <- can.phi * z.star.lag1
  can.sd      <- sqrt(1 - can.phi^2)

  # likelihood impacted by phi.z does not include first day
  R <- sum(dnorm(z.star[, -1], can.mean, can.sd, log=T)) -
       sum(dnorm(z.star[, -1], cur.mean, cur.sd, log=T)) +
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi     <- can.phi
  } }

  results <- list(z=z, acc=acc, att=att, zg=zg,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi)

  return(results)
}