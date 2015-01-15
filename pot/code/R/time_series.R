ts.sample.z <- function(z.star, acc.z, att.z, mh.z, zg,
                        phi, acc.phi, att.phi, mh.phi,
                        y, z.alpha, x.beta, tau, taug, g, prec.cor) {
  nt <- ncol(y)
  nknots <- nrow(z.star)
  for (t in 1:nt) {
  	cur.res   <- y[, t] - x.beta[, t] + z.alpha * zg[, t]
    cur.rss   <- quad.form(prec.cor, sqrt(taug[, t]) * cur.res)
    
    for (k in 1:nknots) {
      att.z[k, t]   <- att.z[k, t] + 1
      can.z.star    <- z.star[, t]
      can.z.star[k] <- rnorm(1, z.star[k, t], mh.z[k, t])  # will be nknots long
      can.z         <- abs(can.z.star)
      can.zg        <- can.z[g[, t]]  # will be ns long
    
      can.res <- y[, t] - x.beta[, t] + z.alpha * can.zg
      can.rss <- quad.form(prec.cor, sqrt(taug[, t]) * can.res)
      
      # prior 
      if (t == 1) {
        mean <- 0
        sd   <- sqrt(1 / tau[k, t])
      } else {
        mean <- phi * z.star[k, (t - 1)]
        sd   <- sqrt((1 - phi^2) / tau[k, t])
      }
    
      R <- -0.5 * sum(can.rss - cur.rss) + 
            dnorm(can.z.star[k], mean, sd, log=T) -
            dnorm(z.star[k, t], mean, sd, log=T)
      
      if (t < nt) {
      	sd.next <- sqrt((1 - phi^2) / tau[k, (t + 1)])
      	R <- R + dnorm(z.star[k, (t + 1)], (phi * can.z.star[k]), sd.next, log=T) -
                 dnorm(z.star[k, (t + 1)], (phi * z.star[k, t]), sd.next, log=T)
      }

      if (!is.na(R)) { if (log(runif(1)) < R) {
        acc.z[k, t]  <- acc.z[k, t] + 1
        z.star[k, t] <- can.z.star[k]
        zg[, t]      <- can.zg
        cur.res      <- can.res
        cur.rss      <- can.rss
      } }
    }
  }
  
  # phi.z
  att.phi  <- att.phi + 1
  z.lag1   <- z.star[, -nt]  # For the mean, we don't need the last day
  cur.mean <- phi * z.lag1
  
  phi.con     <- qnorm(phi)  # transform to R
  can.phi.con <- rnorm(1, phi.con, mh.phi)  # draw candidate
  can.phi     <- pnorm(can.phi.con)  # transform back to (0, 1)
  if (can.phi == 0) {  # numerical stability
  	can.phi <- 0.000001
  } else if (can.phi == 1) {
  	can.phi <- 0.999999
  }
  can.mean <- can.phi * z.lag1  # will be nknots x nt
  R <- sum(dnorm(z.star[, -1], can.mean, sqrt((1 - can.phi^2)/tau), log=T)) -  # tau is inv var
       sum(dnorm(z.star[, -1], cur.mean, sqrt((1 - phi^2)/tau), log=T)) +
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi     <- can.phi
  } }
    
  results <- list(z.star=z.star, acc.z=acc.z, att.z=att.z, zg=zg,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi)
  return(results)
}