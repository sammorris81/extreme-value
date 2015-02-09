updateRhoNu <- function(rho, logrho.m, logrho.s, fixnu, nu, lognu.m, lognu.s,
                        d, gamma, res, taug, prec, cor, logdet.prec, cur.rss,
                        rho.upper=Inf, nu.upper=Inf,
                        att.rho, acc.rho, mh.rho, att.nu, acc.nu, mh.nu) {
  nt <- ncol(res)
  att.rho <- att.rho + 1
  att.nu  <- att.nu + 1

  logrho <- log(rho)
  if (rho.upper == Inf) {
    upper.logrho <- 1
  } else {
    upper.logrho <- pnorm(log(rho.upper), logrho, mh.rho)
  }
  can.rho.u  <- runif(1) * upper.logrho
  can.logrho <- logrho + mh.rho * qnorm(can.rho.u)
  can.rho    <- exp(can.logrho)

  lognu  <- log(nu)
  if (!fixnu) {
    if (nu.upper == Inf) {
      upper.lognu <- 1
    } else {
      upper.lognu <- pnorm(log(nu.upper), lognu, mh.nu)
    }
    can.nu.u  <- runif(1) * upper.lognu
    can.lognu <- lognu + mh.nu * qnorm(can.nu.u)
  } else {
    can.lognu <- lognu
  }
  can.nu <- exp(can.lognu)

  # CorFx gives can.C too, but want a way to pass through gammaless
  # correlation matrix to speed up gamma update
  can.cor     <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, can.rho),
                               error.var=0, smoothness=can.nu, finescale.var=0)
  can.C       <- gamma * can.cor
  diag(can.C) <- 1
  # can.C <- CorFx(d=d, gamma=gamma, rho=can.rho, nu=can.nu)
  can.CC <- tryCatch(chol.inv(can.C, inv=T, logdet=T),
                     error = function(e) {
                       tryCatch(eig.inv(can.C, inv=T, logdet=T, mtx.sqrt=T),
                                error = function(e) {
                                  print(paste("can.rho =", can.rho))
                                  print(paste("can.nu =", can.nu))
                                })
                     })
  can.prec        <- can.CC$prec
  can.logdet.prec <- can.CC$logdet.prec  # this is the sqrt of logdet.prec

  can.rss <- sum(rss(prec=can.prec, y=sqrt(taug) * res))

  R <- -0.5 * (can.rss - cur.rss) +
        nt * (can.logdet.prec - logdet.prec) +
        dnorm(can.logrho, logrho.m, logrho.s, log=T) -
        dnorm(logrho, logrho.m, logrho.s, log=T)

  if (upper.logrho < 1) {  # candidate is not symmetric
    R <- R + dnorm(logrho, can.logrho, mh.rho, log=T) -
             pnorm(upper.logrho, can.logrho, mh.rho, log.p=T) -
             dnorm(can.logrho, logrho, mh.rho, log=T) +
             pnorm(upper.logrho, logrho, mh.rho, log.p=T)
  }

  if (!fixnu) {
    R <- R + dnorm(can.lognu, lognu.m, lognu.s, log=T) -
             dnorm(lognu, lognu.m, lognu.s, log=T)
    if (upper.lognu < 1) {  # candidate is not symmetric
      R <- R + dnorm(lognu, can.lognu, mh.nu, log=T) -
               pnorm(upper.lognu, can.lognu, mh.nu, log.p=T) -
               dnorm(can.lognu, lognu, mh.nu, log=T) +
               pnorm(upper.lognu, lognu, mh.nu, log.p=T)
    }
  }

  if (!is.na(R)) { if (log(runif(1)) < R) {
      rho         <- can.rho
      nu          <- can.nu
      prec        <- can.prec
      cor         <- can.cor
      logdet.prec <- can.logdet.prec
      acc.rho     <- acc.rho + 1
      acc.nu      <- acc.nu + 1
      cur.rss     <- can.rss
  }}

  results <- list(rho=rho, nu=nu, prec=prec, cor=cor, logdet.prec=logdet.prec,
                  cur.rss=cur.rss, att.rho=att.rho, acc.rho=acc.rho,
                  att.nu=att.nu, acc.nu=acc.nu)
}


updateGamma <- function(gamma, gamma.m, gamma.s, d, rho, nu, taug, res, prec,
                        cor, logdet.prec, cur.rss, att, acc, mh) {
  nt <- ncol(res)
  att <- att + 1

  gamma.star     <- transform$probit(gamma)
  can.gamma.star <- rnorm(1, gamma.star, mh)
  can.gamma      <- transform$inv.probit(can.gamma.star)

  can.C  <- can.gamma * cor
  diag(can.C) <- 1
  can.CC <- tryCatch(chol.inv(can.C, inv=T, logdet=T),
                     error = function(e) {
                        tryCatch(eig.inv(can.C, inv=T, logdet=T, mtx.sqrt=T),
                                 error = function(e) {
                                   print(paste("can.gamma =", can.gamma))
                                 })
                     })
  can.prec        <- can.CC$prec
  can.logdet.prec <- can.CC$logdet.prec  # sqrt of logdet.prec

  can.rss <- sum(rss(prec=can.prec, y=sqrt(taug) * res))

  R <- -0.5 * (can.rss - cur.rss) + nt * (can.logdet.prec - logdet.prec) +
        dnorm(can.gamma.star, log=T) - dnorm(gamma.star, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc         <- acc + 1
    gamma       <- can.gamma
    prec        <- can.prec
    logdet.prec <- can.logdet.prec
    cur.rss     <- can.rss
  }}

  results <- list(gamma=gamma, prec=prec, logdet.prec=logdet.prec,
                  cur.rss=cur.rss, acc=acc, att=att)

  return(results)
}

updateRhoNuGamma <- function(rho, logrho.m, logrho.s, fixnu, nu, lognu.m,
                             lognu.s, d, gamma, res, taug, prec,
                             logdet.prec, cur.rss, rho.upper=Inf, nu.upper=Inf,
                             att.rho, acc.rho, mh.rho,
                             att.nu, acc.nu, mh.nu,
                             att.gamma, acc.gamma, mh.gamma) {
  nt <- ncol(res)
  att.rho <- att.rho + 1
  att.nu  <- att.nu + 1
  att.gamma <- att.gamma + 1

  logrho <- log(rho)
  if (rho.upper == Inf) {
    upper.logrho <- 1
  } else {
    upper.logrho <- pnorm(log(rho.upper), logrho, mh.rho)
  }
  can.rho.u  <- runif(1) * upper.logrho
  can.logrho <- logrho + mh.rho * qnorm(can.rho.u)
  can.rho    <- exp(can.logrho)

  lognu  <- log(nu)
  if (!fixnu) {
    if (nu.upper == Inf) {
      upper.lognu <- 1
    } else {
      upper.lognu <- pnorm(log(nu.upper), lognu, mh.nu)
    }
    can.nu.u  <- runif(1) * upper.lognu
    can.lognu <- lognu + mh.nu * qnorm(can.nu.u)
  } else {
    can.lognu <- lognu
  }
  can.nu <- exp(can.lognu)

  gamma.star     <- transform$probit(gamma)
  can.gamma.star <- rnorm(1, gamma.star, mh.gamma)
  can.gamma      <- transform$inv.probit(can.gamma.star)

  # CorFx gives can.C too, but want a way to pass through gammaless
  # correlation matrix to speed up gamma update
  can.C <- CorFx(d=d, gamma=can.gamma, rho=can.rho, nu=can.nu)
  can.CC <- tryCatch(chol.inv(can.C, inv=T, logdet=T),
                     error = function(e) {
                       tryCatch(eig.inv(can.C, inv=T, logdet=T, mtx.sqrt=T),
                                error = function(e) {
                                  print(paste("can.rho =", can.rho))
                                  print(paste("can.nu =", can.nu))
                                })
                     })
  can.prec        <- can.CC$prec
  can.logdet.prec <- can.CC$logdet.prec  # this is the sqrt of logdet.prec

  can.rss <- sum(rss(prec=can.prec, y=sqrt(taug) * res))

  R <- -0.5 * (can.rss - cur.rss) +
        nt * (can.logdet.prec - logdet.prec) +
        dnorm(can.logrho, logrho.m, logrho.s, log=T) -
        dnorm(logrho, logrho.m, logrho.s, log=T) +
        dnorm(can.gamma.star, log=T) - dnorm(gamma.star, log=T)

  if (upper.logrho < 1) {  # candidate is not symmetric
    R <- R + dnorm(logrho, can.logrho, mh.rho, log=T) -
             pnorm(upper.logrho, can.logrho, mh.rho, log.p=T) -
             dnorm(can.logrho, logrho, mh.rho, log=T) +
             pnorm(upper.logrho, logrho, mh.rho, log.p=T)
  }

  if (!fixnu) {
    R <- R + dnorm(can.lognu, lognu.m, lognu.s, log=T) -
             dnorm(lognu, lognu.m, lognu.s, log=T)
    if (upper.lognu < 1) {  # candidate is not symmetric
      R <- R + dnorm(lognu, can.lognu, mh.nu, log=T) -
               pnorm(upper.lognu, can.lognu, mh.nu, log.p=T) -
               dnorm(can.lognu, lognu, mh.nu, log=T) +
               pnorm(upper.lognu, lognu, mh.nu, log.p=T)
    }
  }

  if (!is.na(R)) { if (log(runif(1)) < R) {
      rho         <- can.rho
      nu          <- can.nu
      gamma       <- can.gamma
      prec        <- can.prec
      logdet.prec <- can.logdet.prec
      acc.rho     <- acc.rho + 1
      acc.nu      <- acc.nu + 1
      acc.gamma   <- acc.gamma + 1
      cur.rss     <- can.rss
  }}

  results <- list(rho=rho, nu=nu, gamma=gamma, prec=prec,
                  logdet.prec=logdet.prec, cur.rss=cur.rss, att.rho=att.rho,
                  acc.rho=acc.rho, att.nu=att.nu, acc.nu=acc.nu,
                  att.gamma=att.gamma, acc.gamma=acc.gamma)
}