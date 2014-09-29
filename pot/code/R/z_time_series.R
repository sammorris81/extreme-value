# zstar is nknots x nt - actual values for time series
ts.sample.z <- function(zstar, acc.z, att.z, mh.z, 
                        phi, acc.phi, att.phi, mh.phi,
                        y, z.alpha, x.beta, tau, taug, g, prec.cor) {
  nt <- ncol(y)
  
  for (t in 1:nt) {
  	att.z[, t] <- att.z[, t] + 1
    can.zstar <- rnorm(nknots, zstar[, t], mh.z[, t])  # will be nknots long
    can.z <- abs(can.zstar)
    can.zg <- can.z[g[, t]]  # will be ns long
    can.res <- y - x.beta[, t] + z.alpha * can.zg
    cur.res <- y - x.beta[, t] + z.alpha * abs(zstar[g[, t], t])
    can.rss <- quad.form(prec.cor, sqrt(taug[, t]) * can.res)
    cur.rss <- quad.form(prec.cor, sqrt(taug[, t]) * cur.res)
    
    # prior 
    if (t == 1) {
      mean <- 0
      sd <- 1
    } else {
      mean <- phi * zstar[, (t - 1)]
      sd <- sqrt((1 - phi^2) / tau[, t])
    }
    
    R <- -0.5 * sum(can.rss - cur.rss) + 
          sum(dnorm(can.zstar, mean, sd, log=TRUE)) -
          sum(dnorm(zstar[, t], mean, sd, log=TRUE))
    
    if (!is.na(R)) { if (log(runif(1)) < R) {
      acc.z[, t] <- acc.z[, t] + 1
      zstar[, t] <- can.zstar
    } }
    
  }
  
  # phi.z
  att.phi <- att.phi + 1
  z.lag1 <- cbind(rep(0, nknots), zstar[, -(nt + 1), drop=F])
  cur.mean <- phi * z.lag1
  
  phi_con <- qnorm((phi + 1) / 2)  # transform to R
  can_phi_con <- rnorm(1, phi_con, mh.phi)  # draw candidate
  can.phi <- 2 * pnorm(can_phi_con) - 1  # transform back to (-1, 1)
  can.mean <- can.phi * z.lag1  # will be nknots x nt
  R <- sum(dnorm(zstar, can.mean, sqrt(1 - can.phi^2)/tau, log=T)) - 
       sum(dnorm(zstar, cur.mean, sqrt(1 - phi^2)/tau, log=T)) +
       dnorm(phi_con, log=T) - dnorm(phi, log=T)
  
  if (!is.nar(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi <- can.phi
  } }
  
  results <- list(zstar=zstar, phi=phi)
  return(results)
}

