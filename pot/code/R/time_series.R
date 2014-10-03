ts.sample.z <- function(z, acc.z, att.z, mh.z, zg,
                        phi, acc.phi, att.phi, mh.phi,
                        y, z.alpha, x.beta, tau, taug, g, prec.cor) {
  nt <- ncol(y)
  nknots <- nrow(z)
  logz <- log(z)
  for (t in 1:nt) {
  	att.z[, t] <- att.z[, t] + 1
    can.logz <- rnorm(nknots, logz[, t], mh.z[, t])  # will be nknots long
    can.z <- exp(can.logz)
    can.zg <- can.z[g[, t]]  # will be ns long
    can.res <- y - x.beta[, t] + z.alpha * can.zg
    cur.res <- y - x.beta[, t] + z.alpha * zg[, t]
    can.rss <- quad.form(prec.cor, sqrt(taug[, t]) * can.res)
    cur.rss <- quad.form(prec.cor, sqrt(taug[, t]) * cur.res)

    # prior 
    if (t == 1) {
      mean <- 0
      sd <- sqrt(1 / tau[, t])
    } else {
      mean <- phi * z[, (t - 1)]
      sd <- sqrt((1 - phi^2) / tau[, t])
    }
    
    R <- -0.5 * sum(can.rss - cur.rss) + 
          sum(log(dfoldnorm(can.z, mean, sd))) -
          sum(log(dfoldnorm(z[, t], mean, sd)))

    if (!is.na(R)) { if (log(runif(1)) < R) {
      acc.z[, t] <- acc.z[, t] + 1
      z[, t]  <- can.z
      zg[, t] <- can.zg
    } }
    
  }
  
  # phi.z
  att.phi <- att.phi + 1
  z.lag1 <- cbind(rep(0, nknots), z[, -nt, drop=F])
  cur.mean <- phi * z.lag1
  
  phi.con     <- qnorm(phi)  # transform to R
  can.phi.con <- rnorm(1, phi.con, mh.phi)  # draw candidate
  can.phi     <- pnorm(can.phi.con)  # transform back to (0, 1)
  if (can.phi == 0) {  # numerical stability
  	can.phi = 0.000001
  } else if (can.phi == 1) {
  	can.phi = 0.999999
  }
  can.mean    <- can.phi * z.lag1  # will be nknots x nt
  R <- sum(log(dfoldnorm(z, can.mean, sqrt((1 - can.phi^2)/tau)))) -  # tau is inv var
       sum(log(dfoldnorm(z, cur.mean, sqrt((1 - phi^2)/tau)))) +
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi <- can.phi
  } }
  
  results <- list(z=z, acc.z=acc.z, att.z=att.z, zg=zg,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi)
  return(results)
}

ts.sample.tau <- function(tau, acc.tau, att.tau, mh.tau, taug,
                          phi, att.phi, acc.phi, mh.phi,
                          s, s.a, s.b, 
                          res, prec.cor, g) {
  
  nt <- ncol(tau)
  nknots <- nrow(tau)
  # update tau terms
  for (t in 1:nt) {
    res.t <- res[, t]
    curll.t <- 0.5 * sum(log(taug[, t])) - 
               0.5 * quad.form(prec.cor, sqrt(taug[, t]) * res.t)
    
    for (k in 1:nknots) {
      these <- which(g[, t] == k)
      nparts <- length(these)
      
      if (t == 1) {
      	mean <- 0
      	sd <- s
      } else {
      	mean <- phi * log(tau[k, (t - 1)])
      	sd <- s * sqrt(1 - phi^2)
      }
      
      if (nparts == 0) {  # we just draw from the prior because no likelihood
        logtau <- rnorm(1, mean, sd)
        tau[k, t] <- exp(logtau)
      } else {  # do a MH update
      	att.tau[k, t] <- att.tau[k, t] + 1
        canlogtau <- log(tau[, t])  # pull out all taus for a day
        canlogtau[k] <- rnorm(1, log(tau[k, t]), mh.tau)  # get candidate for knot k
        canlogtaug <- canlogtau[g[, t]]  # transform to length ns
      	
        canll <- 0.5 * sum(canlogtau) - 
                 0.5 * quad.form(prec.cor, sqrt(exp(canlogtaug)) * res.t)
                 
        R <- canll - curll.t +
             dnorm(canlogtau[k], mean, sd, log=TRUE) -
             dnorm(log(tau[k, t]), mean, sd, log=TRUE)
             
        if (!is.na(R)) { if (log(runif(1)) < R) {
          acc.tau[k, t] <- acc.tau[k, t] + 1
          tau[, t] <- exp(canlogtau)
          taug[, t] <- exp(canlogtaug)
          curll.t <- canll
        } }
      }  # fi nparts == 0
    }  # end k in 1:nknots
  }  # end t in 1:nt
  
  # phi.tau
  att.phi <- att.phi + 1
  logtau.lag1 <- cbind(rep(0, nknots), log(tau[, -nt, drop=F]))  # should be nknots x nt
  cur.mean <- phi * logtau.lag1
  
  phi.con <- qnorm((phi + 1) / 2)  # transform to R
  can.phi.con <- rnorm(1, phi.con, mh.phi)  # draw candidate
  can.phi <- 2 * pnorm(can.phi.con) - 1  # transform back to (-1, 1)
  can.mean <- can.phi * logtau.lag1
  
  R <- sum(dnorm(log(tau), can.mean, sqrt(1 - can.phi^2) * s, log=T)) - 
       sum(dnorm(log(tau), cur.mean, sqrt(1 - can.phi^2) * s, log=T)) + 
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)
       
  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi <- can.phi
  } }
  
  # s: prior is IG(a, b)
  a.star <- s.a + nt * nknots / 2
  b.star <- s.b
  for (t in 1:nt) {
    if (t == 1) {
      b.star <- b.star + sum(log(tau[, t])^2 / 2)
    } else {
      b.star <- b.star + sum((log(tau[, t]) - phi * log(tau[, (t - 1)]))^2 / (2 * (1 - phi^2)))
    }
  }
  s <- 1 / rgamma(1, a.star, b.star)
  
  results <- list(tau=tau, att.tau=att.tau, acc.tau=acc.tau,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi,
                  s)
  return(results)                         
}
