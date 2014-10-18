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
  z.lag1   <- cbind(rep(0, nknots), z.star[, -nt, drop=F])
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
  R <- sum(dnorm(z.star, can.mean, sqrt((1 - can.phi^2)/tau), log=T)) -  # tau is inv var
       sum(dnorm(z.star, cur.mean, sqrt((1 - phi^2)/tau), log=T)) +
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi     <- can.phi
  } }
    
  results <- list(z.star=z.star, acc.z=acc.z, att.z=att.z, zg=zg,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi)
  return(results)
}

ts.sample.tau <- function(tau, acc.tau, att.tau, mh.tau, 
                          att.tau.ns, acc.tau.ns, mh.tau.ns, taug,
                          phi, att.phi, acc.phi, mh.phi,
                          s, s.a, s.b, tau.alpha, tau.beta,
                          res, prec.cor, g, z) {
  
  nt       <- ncol(tau)
  nknots   <- nrow(tau)
  tau.star <- cop.inv.IG(tau=tau, phi=phi, alpha=tau.alpha, beta=tau.beta)
  # update tau terms
  for (t in 1:nt) {
    res.t <- res[, t]
    # get the current day's likelihood
    curll <- 0.5 * sum(log(taug[, t])) - 
             0.5 * quad.form(prec.cor, sqrt(taug[, t]) * res.t) + 
             0.5 * sum(log(tau[, t])) - 
             0.5 * sum(tau[, t] * (z[, t])^2)
    
    for (k in 1:nknots) {
      these  <- which(g[, t] == k)
      nparts <- length(these)
      
      if (t == 1) {
      	mean <- 0
      	sd   <- sqrt(s)
      } else {
      	mean <- phi * tau.star[k, (t - 1)]
      	sd   <- sqrt(s * (1 - phi^2))
      }
      
      # draw candidate tau.star from Normal 
      # prior is on tau.stars
      # likelihood uses tau
      
      if (nparts == 0) {  # we just draw from the prior because no likelihood
        tau.star[k, t] <- rnorm(1, mean, sd)
        res.std        <- (tau.star[k, t] - mean) / sd
        tau[k, t]      <- qgamma(pnorm(res.std), shape=tau.alpha, rate=tau.beta)
      } else {  # do a MH update - First try with block by day
        # candidate moves are in the AR(1) distribution
        att.tau.ns[nparts] <- att.tau.ns[nparts] + 1
        att.tau[k, t]      <- att.tau[k, t] + 1
      	can.tau.star       <- tau.star[, t]  # pull out all taus for a day
        can.tau.star[k]    <- rnorm(1, tau.star[k, t], mh.tau.ns[nparts])  # get candidate for knot k
        # can.tau.star[k] <- rnorm(1, tau.star[k, t], mh.tau[k, t])  # get candidate for knot k
        can.taug.star      <- can.tau.star[g[, t]]  # transform to length ns
        
        # transform to IG marginals
        can.tau    <- tau[, t]
        res.std    <- (can.tau.star[k] - mean)/sd
        can.tau[k] <- qgamma(pnorm(res.std), shape=tau.alpha, rate=tau.beta)
        can.taug   <- can.tau[g[, t]]
      	
        canll <- 0.5 * sum(log(can.taug)) - 
                 0.5 * quad.form(prec.cor, sqrt(can.taug) * res.t) +
                 0.5 * sum(log(can.tau)) - 
                 0.5 * sum(can.tau * (z[, t])^2)
               
        R <- canll - curll +
             dnorm(can.tau.star[k], mean, sd, log=TRUE) -
             dnorm(tau.star[k, t], mean, sd, log=TRUE)
        
        if (t < nt) {
          sd.next <- sqrt(s * (1 - phi^2))
          R <- R + dnorm(tau.star[k, (t + 1)], (phi * can.tau.star[k]), sd.next, log=T) -
                   dnorm(tau.star[k, (t + 1)], (phi * tau.star[k, t]), sd.next, log=T)
        }
             
        if (!is.na(R)) { if (log(runif(1)) < R) {
          acc.tau.ns[nparts] <- acc.tau.ns[nparts] + 1
          acc.tau[k, t]      <- acc.tau[k, t] + 1
          tau[k, t]          <- can.tau[k]
          tau.star[k, t]     <- can.tau.star[k]
          taug[these, t]     <- can.taug[these]
          curll              <- canll
        } }
       
      }  # fi nparts == 0
    }  # end k in 1:nknots
  }  # end t in 1:nt
  
  # phi.tau
  att.phi     <- att.phi + 1
  tau.star.lag1 <- cbind(rep(0, nknots), tau.star[, -nt, drop=F])  # should be nknots x nt
  cur.mean    <- phi * tau.star.lag1
  
  phi.con     <- qnorm((phi + 1) / 2)  # transform to R
  can.phi.con <- rnorm(1, phi.con, mh.phi)  # draw candidate
  can.phi     <- 2 * pnorm(can.phi.con) - 1  # transform back to (-1, 1)
  can.mean    <- can.phi * tau.star.lag1
  
  R <- sum(dnorm(tau.star, can.mean, sqrt(s * (1 - can.phi^2)), log=T)) - 
       sum(dnorm(tau.star, cur.mean, sqrt(s * (1 - can.phi^2)), log=T)) + 
       dnorm(can.phi.con, log=T) - dnorm(phi.con, log=T)
       
  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc.phi <- acc.phi + 1
    phi     <- can.phi
  } }
  
  # s: prior is IG(a, b)
  a.star <- s.a + nt * nknots / 2
  b.star <- s.b
  for (t in 1:nt) {
    if (t == 1) {
      b.star <- b.star + sum(tau.star[, t]^2 / 2)
    } else {
      b.star <- b.star + sum(tau.star[, t] - 
                             phi * tau.star[, (t - 1)]^2 / (2 * (1 - phi^2)))
    }
  }
  s <- 1 / rgamma(1, a.star, b.star)
  
  results <- list(tau=tau, att.tau=att.tau, acc.tau=acc.tau,
                  phi=phi, att.phi=att.phi, acc.phi=acc.phi,
                  s)
  return(results)                         
}
