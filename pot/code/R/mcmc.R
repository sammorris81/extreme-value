################################################################################
# MCMC
#
# TODO: Add in model description here
#
#ASSUMES x and y are in [0,1]^2
################################################################################
if (!exists("conditional.Rcpp")) {
  source('condmean_cpp.R')
}
if (!exists("g.Rcpp")) {
  source('mem_cpp.R')
}
if (!exists("z.Rcpp")) {
  source('z_update_cpp.R')
}
source('imputeY.R')
source('updateBeta.R')
source('updateCov.R')
source('updateKnots.R')
source('updateLambda.R')
source('updatePhi.R')
source('updateTau.R')
source('updateZ.R')
source('predictY.R')

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL,
                 min.s, max.s,  # don't want to specify defaults
                 thresh.all=0, thresh.quant=T, nknots=1, keep.knots=F,
                 iters=5000, burn=1000, update=100, thin=1,
                 iterplot=F, plotname=NULL, method="t",
                 # just to debug temporal parts. eventually change to temporal=F
                 temporalw=F, temporaltau=F, temporalz=F,
                 # initial values
                 beta.init=NULL, tau.init=1,
                 tau.alpha.init=0.1, tau.beta.init=0.1,
                 rho.init=5, nu.init=0.5, gamma.init=0.5,
                 # priors
                 beta.m=0, beta.s=10,
                 tau.alpha.m=0, tau.alpha.s=1,
                 tau.beta.a=1, tau.beta.b=1,
                 logrho.m=0, logrho.s=1, rho.upper=NULL,
                 lognu.m=-1.2, lognu.s=1, nu.upper=NULL,
                 gamma.m=0, gamma.s=1,
                 # covariance model
                 cov.model="matern",  # or "exponential"
                 rho.prior="cont",  # or "disc"
                 # skew inits
                 z.init=1, lambda.init=0.5,
                 # skew priors
                 lambda.a=1, lambda.b=1, skew=T,
                 thresh.site.specific=F, thresh.site=NULL,
                 # troubleshooting
                 debug=F, fixhyper=F, tau.t, z.t, fixknots=F, knots.init=NULL,
                 fixbeta=F
        ){

  library(SpatialTools)
  library(fields)
  library(emulator)
  start.time <- proc.time()

  ##############################################
  # Initial setup
  ##############################################
  ns <- nrow(y)    # number of sites
  nt <- ncol(y)    # number of days
  p  <- dim(x)[3]  # number of covariates

  predictions <- !is.null(s.pred) & !is.null(x.pred)
  np <- 0
  y.pred <- NULL
  if (predictions) {
    np <- nrow(s.pred)
    d12 <- rdist(s.pred, s)
    d11 <- rdist(s.pred, s.pred)
    diag(d11) <- 0
    y.pred <- array(0, c((iters-burn), np, nt))
  }

  d       <- rdist(s)  # distance between sites
  diag(d) <- 0

  y.init <- y  # store the initial y values

  # store the loc/day for observed values below thresh.
  if (thresh.quant & thresh.all >= 0) {  # threshold based on sample quantiles
    if (thresh.all < 0 | thresh.all > 1) {
      warning("Warning: quantile for thresholding must be between 0 and 1.")
    }
    thresh.all.q <- quantile(y, thresh.all, na.rm=T)
  } else {
    thresh.all.q <- thresh.all
  }
  thresh.all.mtx <- matrix(thresh.all.q, ns, nt)  # matrix for easy replacement
  if (thresh.site.specific) {
    if (is.null(thresh.site)) {
      warning("Warning: setting site-specific time series threshold to thresh.")
      thresh.site <- thresh.all
    }
    if ((thresh.site < 0) | (thresh.site > 1)) {
      stop("Error: thresh.site should be the desired site-specific quantile
            between 0 and 1")
    }
    # if it's site specific, then we want to keep everything over the 95th
    # quantile for the data set and the max for each site in a matrix that's
    # ns x nt
    thresh.site <- apply(y, 1, quantile, probs=thresh.site, na.rm=T)
    thresh.site.mtx <- matrix(thresh.site, ns, nt)
    thresh.mtx <- ifelse(
      thresh.site < thresh.all.mtx,
      thresh.site,
      thresh.all.mtx
    )
  } else {  # if the mean doesn't have a site component, same threshold for all
    cat("\t no site-specific threshold set \n")
    thresh.mtx  <- thresh.all.mtx
  }
  thresh.obs  <- !is.na(y) & (y < thresh.mtx)
  if (sum(thresh.obs) > 0) {
    thresholded <- T
  } else {
    thresholded <- F
  }
  # y[thresh.obs] <- thresh.mtx[thresh.obs] / 2

  missing.obs <- is.na(y)
  if (sum(missing.obs) > 0) {
    y[missing.obs]  <- mean(y, na.rm=T)
    missing <- T
  } else {
    missing <- F
  }

  # initialize partition
  if (!fixknots) {
    if (is.null(knots.init)) {
      # constrain the initial draw to be near the middle for stability
      lower.1 <- min.s[1] + (max.s[1] - min.s[1]) / 3
      upper.1 <- min.s[1] + 2 * (max.s[1] - min.s[1]) / 3
      lower.2 <- min.s[2] + (max.s[2] - min.s[2]) / 3
      upper.2 <- min.s[2] + 2 * (max.s[2] - min.s[2]) / 3
      knots <- array(NA, dim=c(nknots, 2, nt))
      knots[, 1, ] <- runif(nknots * nt, lower.1, upper.1)
      knots[, 2, ] <- runif(nknots * nt, lower.2, upper.2)
    } else {
      knots <- array(knots.init, dim=c(nknots, 2, nt))
    }
    knots.star        <- array(NA, c(nknots, 2, nt))
    knots.star[, 1, ] <- transform$probit(knots[, 1, ], lower=min.s[1],
                                          upper=max.s[1])
    knots.star[, 2, ] <- transform$probit(knots[, 2, ], lower=min.s[2],
                                          upper=max.s[2])
  } else {
    if (is.null(knots.init)) {
      stop("You must specify knots.init to fix the knots")
    }
    knots <- knots.init
  }

  # initialize parameters
  if (is.null(beta.init)) {
    beta    <- rep(0, p)
    # beta[1] <- mean(y)
  } else {
    beta <- beta.init
  }
  x.beta  <- matrix(beta[1], ns, nt)

  # initialize partitioning
  if (nknots == 1) {
    g <- matrix(1, ns, nt)
  } else {
    g <- matrix(0, ns, nt)
    for (t in 1:nt) {
      g[, t] <- mem(s, knots[, , t])
    }
  }

  # initialize variance
  taug <- matrix(0, ns, nt)
  if (method == "gaussian") {  # single knot for all days
  	if (length(tau.init) > 1) {
  	  stop("for gaussian, tau.init should be a single value")
  	}
  	tau <- matrix(tau.init, nknots, nt)
    taug <- matrix(tau.init, ns, nt)
  } else if (method == "t") {  # knots vary by day and partition
  	if (length(tau.init) == 1) {
  	  cat("\t initializing all tau terms to", tau.init, "\n")
  	}
  	tau  <- matrix(tau.init, nknots, nt)
  	for (t in 1:nt) {
      taug[, t] <- tau[g[, t], t]
    }
  }
  if (debug) {
    tau <- matrix(tau.t, nknots, nt)
    for (t in 1:nt) {
  	  taug[, t] <- tau[g[, t], t]
    }
  }

  zg <- matrix(0, ns, nt)
  if (length(z.init) == 1 && skew) {
    cat("\t initializing all z terms to", z.init, "\n")
  }
  z <- matrix(z.init, nknots, nt)
  for (t in 1:nt) {
    zg[, t] <- z[g[, t], t]
  }
  if (debug) {
    z <- matrix(z.t, nknots, nt)
    for (t in 1:nt) {
      zg[, t] <- z[g[, t], t]
    }
  }

  if (skew) {
    lambda <- lambda.init
    if (lambda < 0) {
      lambda.1 <- -1
    } else {
      lambda.1 <- 1
    }
    lambda.2 <- 1 / (lambda)^2
  } else {
    lambda <- lambda.1 <- 0
    lambda.2 <- 0
  }

  # easier to keep calculations in the precision scale for MCMC
  sigma2    <- 1 / tau
  sigma2g   <- 1 / taug
  tau.alpha <- tau.alpha.init
  tau.beta  <- tau.beta.init

  # initialize spatial covariance
  rho    <- rho.init
  logrho <- log(rho)
  if (is.null(rho.upper)) {
    rho.upper <- Inf
  }
  nu     <- nu.init
  lognu  <- log(nu)
  if (is.null(nu.upper)) {
    nu.upper <- Inf
  }
  gamma  <- gamma.init

  fixnu <- F
  if (cov.model == "exponential") {
    nu <- 0.5
    fixnu <- T
  }

  #### Precomputed initialized values.
  # setup correlation matrix without multiplying by gamma
  # simple.cov.sp only needs to be updated for rho and nu updates
  cor <- simple.cov.sp(D=d, sp.type="matern", sp.par=c(1, rho),
                       error.var=0, smoothness=nu, finescale.var=0)
  C <- gamma * cor
  diag(C) <- 1
  CC <- tryCatch(chol.inv(C, inv=T, logdet=T),
                 error = function(e) {
                   eig.inv(C, inv=T, logdet=T, mtx.sqrt=T)
                 })
  prec <- CC$prec
  logdet.prec <- CC$logdet.prec  # already includes 0.5

  # time series in the random knots/partitions
  if (temporalw) {
    phi.w <- 0
    acc.phi.w <- att.phi.w <- mh.phi.w <- 1
  }
  if (temporaltau) {
    phi.tau <- 0
    acc.phi.tau <- att.phi.tau <- mh.phi.tau <- 1
  } else {
    acc.tau.low <- acc.tau.high <- 0
    tau.trials <- nknots * nt * 3
  }

  acc.tau <- att.tau <- matrix(1, nrow=nknots, ncol=nt)
  mh.tau <- matrix(0.05, nknots, nt)
  nparts.tau <- matrix(1, nrow=nknots, ncol=nt)

  if (temporalz) {
    phi.z <- 0
    acc.z <- att.z <- mh.z <- matrix(1, nknots, nt)
    acc.phi.z <- att.phi.z <- mh.phi.z <- 1
  }

  # MH tuning params
  acc.w      <- att.w      <- mh.w     <- matrix(0.15, nknots, nt)  # knot locs
  acc.delta  <- att.delta  <- mh.delta <- 0.1
  acc.rho    <- att.rho    <- mh.rho   <- 0.1
  acc.nu     <- att.nu     <- mh.nu    <- 0.1
  acc.gamma  <- att.gamma  <- mh.gamma <- 1

  # storage
  keepers.tau       <- array(NA, dim=c(iters, nknots, nt))
  keepers.beta      <- matrix(NA, nrow=iters, ncol=p)
  keepers.ll        <- matrix(NA, nrow=iters, ncol=nt)
  keepers.tau.alpha <- rep(NA, iters)
  keepers.tau.beta  <- rep(NA, iters)
  keepers.nparts    <- matrix(NA, nrow=iters, ncol=nt)
  keepers.rho       <- rep(NA, iters)
  keepers.nu        <- rep(NA, iters)
  keepers.gamma     <- rep(NA, iters)
  keepers.z         <- array(NA, dim=c(iters, nknots, nt))
  keepers.lambda    <- rep(NA, iters)
  keepers.lambda.1 <- keepers.lambda.2 <- rep(NA, iters)
  keepers.avgparts  <- matrix(NA, nrow=iters, ncol=nt)  # avg partitions per day
  if (keep.knots & (nknots > 1)) {
    keepers.knots     <- array(NA, dim=c(iters, nknots, 2, nt))
  }
  if (temporalz) {
    keepers.phi.z <- rep(NA, iters)
  }
  if (temporalw) {
    keepers.phi.w <- rep(NA, iters)
  }
  if (temporaltau) {
    keepers.phi.tau <- rep(NA, iters)
  }
  return.iters      <- (burn + 1):iters

  tic <- proc.time()
  for (iter in 1:iters) { for (ttt in 1:thin) {


    # data imputation
    if (thresholded) {  # do data imputation and store as y
      mu <- x.beta + lambda.1 * zg
      y <- imputeY(y=y, taug=taug, mu=mu, obs=thresh.obs, cor=cor,
                   gamma=gamma, thresh.mtx=thresh.mtx)
    }

    # missing values
    if (missing) {
      mu <- x.beta + lambda.1 * zg
      y <- imputeY(y=y, taug=taug, mu=mu, obs=missing.obs, cor=cor,
                   gamma=gamma)
    }

    # update beta
    beta <- updateBeta(beta.m=beta.m, beta.s=beta.s, x=x, y=y, zg=zg,
                       lambda.1=lambda.1, taug=taug, prec=prec)
    for (t in 1:nt) {
      x.beta[, t] <- x[, t, ] %*% beta
    }

    # update partitions
    if ((nknots > 1) & (!fixknots)) {
      avgparts <- rep(0, nt)
      knots.update <- updateKnotsTS(phi=phi.w, knots=knots, g=g, ts=temporalw,
                                    tau=tau, z=z, s=s, min.s=min.s, max.s=max.s,
                                    x.beta=x.beta, lambda.1=lambda.1, y=y,
                                    prec=prec, att=att.w, acc=acc.w,
                                    mh=mh.w, att.phi=att.phi.w,
                                    acc.phi=acc.phi.w, mh.phi=mh.phi.w)
      knots.star <- knots.update$knots.star
      knots      <- knots.update$knots
      g          <- knots.update$g
      taug       <- knots.update$taug
      zg         <- knots.update$zg
      acc.w      <- knots.update$acc
      att.w      <- knots.update$att
      phi.w      <- knots.update$phi
      acc.phi.w  <- knots.update$acc.phi
      att.phi.w  <- knots.update$att.phi

      # update metropolis sds
      if (iter < (burn / 2)) {
        mh.update  <- mhupdate(acc=acc.w, att=att.w, mh=mh.w,
                               nattempts=nknots * nt)
        acc.w      <- mh.update$acc
        att.w      <- mh.update$att
        mh.w       <- mh.update$mh

        if (temporalw) {
          mh.update <- mhupdate(acc=acc.phi.w, att=att.phi.w, mh=mh.phi.w)
          acc.phi.w <- mh.update$acc
          att.phi.w <- mh.update$att
          mh.phi.w  <- mh.update$mh
        }
      }
    }  # fi nknots > 1

    # update tau
    # all taus require mu and res
    mu <- x.beta + lambda.1 * zg
    res <- y - mu
    if (method == "gaussian") {  # single random effect for all days
      tau.update <- updateTauGaus(res=res, prec=prec,
                                  tau.alpha=tau.alpha, tau.beta=tau.beta)
      tau  <- matrix(tau.update, nknots, nt)
      taug <- matrix(tau.update, ns, nt)

    } else if (method == "t") {
      if (!temporaltau) {
        tau.update <- updateTau(tau=tau, taug=taug, y=y, mu=mu, g=g, res=res,
                                nparts.tau=nparts.tau, prec=prec,
                                z=z, lambda.2=lambda.2,
                                tau.alpha=tau.alpha, tau.beta=tau.beta,
                                skew=skew, obs=thresh.obs,
                                thresh.mtx=thresh.mtx,
                                att=att.tau, acc=acc.tau, mh=mh.tau)
        tau          <- tau.update$tau
        taug         <- tau.update$taug
        acc.tau      <- tau.update$acc
        att.tau      <- tau.update$att

        # # Note: this doesn't really tune the algorithm because mh tau isn't used
        # if (iter < (burn / 2)) {
        #   mh.update  <- mhupdate(acc=acc.tau, att=att.tau, mh=mh.tau)
        #   acc.tau <- mh.update$acc
        #   att.tau <- mh.update$att
        #   mh.tau  <- mh.update$mh
        # }
      } else {
        tau.update <- updateTauTS(phi=phi.tau, tau=tau, taug=taug,
                                  g=g, res=res, nparts.tau=nparts.tau,
                                  prec=prec, z=z, lambda.2=lambda.2,
                                  tau.alpha=tau.alpha,
                                  tau.beta=tau.beta, skew=skew,
                                  att=att.tau, acc=acc.tau, mh=mh.tau,
                                  att.phi=att.phi.tau, acc.phi=acc.phi.tau,
                                  mh.phi=mh.phi.tau)

          tau         <- tau.update$tau
          taug        <- tau.update$taug
          phi.tau     <- tau.update$phi
          acc.tau     <- tau.update$acc
          att.tau     <- tau.update$att
          acc.phi.tau <- tau.update$acc.phi
          att.phi.tau <- tau.update$att.phi


          # update metropolis sds
          if (iter < (burn / 2)) {
            mh.update <- mhupdate(acc=acc.tau, att=att.tau, mh=mh.tau)
            acc.tau   <- mh.update$acc
            att.tau   <- mh.update$att
            mh.tau    <- mh.update$mh

            mh.update  <- mhupdate(acc=acc.phi.tau, att=att.phi.tau,
                                    mh=mh.phi.tau)
            acc.phi.tau <- mh.update$acc
            att.phi.tau <- mh.update$att
            mh.phi.tau  <- mh.update$mh
          }
      }

      # update tau.alpha and tau.beta
      tau.alpha <- updateTauAlpha(tau=tau, tau.beta=tau.beta)
      tau.beta  <- updateTauBeta(tau=tau, tau.alpha=tau.alpha,
                                 tau.beta.a=tau.beta.a, tau.beta.b=tau.beta.b)
      if (fixhyper) {
        tau.alpha <- 3
        tau.beta  <- 8
      }
    }  # fi method == t

    # update rho and nu and gamma
    # update rss with new residuals and taug
    mu <- x.beta + lambda.1 * zg
    res <- y - mu
    cur.rss <- sum(rss(prec=prec, y=sqrt(taug) * res))

    # rho and nu
    rhonu.update <- updateRhoNu(rho=rho, logrho.m=logrho.m, logrho.s=logrho.s,
                                fixnu=fixnu, nu=nu, lognu.m=lognu.m,
                                lognu.s=lognu.s, d=d, gamma=gamma,
                                res=res, taug=taug, prec=prec, cor=cor,
                                logdet.prec=logdet.prec, cur.rss=cur.rss,
                                rho.upper=rho.upper, nu.upper=nu.upper,
                                y=y, mu=mu, obs=thresh.obs,
                                thresh.mtx=thresh.mtx,
                                att.rho=att.rho, acc.rho=acc.rho, mh.rho=mh.rho,
                                att.nu=att.nu, acc.nu=acc.nu, mh.nu=mh.nu)
    rho         <- rhonu.update$rho
    nu          <- rhonu.update$nu
    prec        <- rhonu.update$prec
    cor         <- rhonu.update$cor
    logdet.prec <- rhonu.update$logdet.prec
    cur.rss     <- rhonu.update$cur.rss
    att.rho     <- rhonu.update$att.rho
    acc.rho     <- rhonu.update$acc.rho
    att.nu      <- rhonu.update$att.nu
    acc.nu      <- rhonu.update$acc.nu

    # gamma - make sure the cor argument does not include gamma
    gamma.update <- updateGamma(gamma=gamma, gamma.m=gamma.m, gamma.s=gamma.s,
                                d=d, rho=rho, nu=nu, taug=taug, res=res,
                                prec=prec, y=y, mu=mu, obs=thresh.obs,
                                thresh.mtx=thresh.mtx,
                                cor=cor, logdet.prec=logdet.prec,
                                cur.rss=cur.rss, att=att.gamma, acc=acc.gamma,
                                mh=mh.gamma)
    gamma       <- gamma.update$gamma
    prec        <- gamma.update$prec
    logdet.prec <- gamma.update$logdet.prec
    cur.rss     <- gamma.update$cur.rss
    acc.gamma   <- gamma.update$acc
    att.gamma   <- gamma.update$att

    if (iter < (burn / 2)) {
      mh.update <- mhupdate(acc=acc.rho, att=att.rho, mh=mh.rho)
      acc.rho   <- mh.update$acc
      att.rho   <- mh.update$att
      mh.rho    <- mh.update$mh

      mh.update <- mhupdate(acc=acc.nu, att=att.nu, mh=mh.nu)
      acc.nu    <- mh.update$acc
      att.nu    <- mh.update$att
      mh.nu     <- mh.update$mh

      mh.update <- mhupdate(acc=acc.gamma, att=att.gamma, mh=mh.gamma)
      acc.gamma <- mh.update$acc
      att.gamma <- mh.update$att
      mh.gamma  <- mh.update$mh
    }

    if (fixhyper) {
      nu    <- 0.5
      rho   <- 1
      gamma <- 0.9
    }

    # update skew parameters: lambda and z
    if (skew) {
      lambda.1 <- updateLambda1(x.beta=x.beta, zg=zg, y=y, prec=prec,
                                taug=taug)
      lambda.1 <- 1
      lambda.2 <- updateLambda2(lambda.a=lambda.a, lambda.b=lambda.b,
                                z=z, tau=tau)
      lambda <- lambda.1 / sqrt(lambda.2)

      if (!temporalz) {
        mu <- x.beta + lambda.1 * zg
        z.update <- updateZ(y=y, x.beta=x.beta, zg=zg, prec=prec,
                            tau=tau, mu=mu, taug=taug, g=g,
                            lambda.1=lambda.1, lambda.2=lambda.2)

        z  <- z.update$z
        zg <- z.update$zg
      } else {
        z.update <- updateZTS(z=z, zg=zg, y=y, lambda.1=lambda.1,
                              lambda.2=lambda.2, x.beta=x.beta,
                              phi=phi.z, tau=tau, taug=taug, g=g,
                              prec=prec, acc=acc.z, att=att.z, mh=mh.z,
                              acc.phi=acc.phi.z, att.phi=att.phi.z,
                              mh.phi=mh.phi.z)
        z      <- z.update$z
        zg     <- z.update$zg
        att.z  <- z.update$att
        acc.z  <- z.update$acc
        phi.z  <- z.update$phi
        att.phi.z <- z.update$att.phi
        acc.phi.z <- z.update$acc.phi

        if (iter < (burn / 2)) {
          mh.update <- mhupdate(acc=acc.z, att=att.z, mh=mh.z)
          acc.z     <- mh.update$acc
          att.z     <- mh.update$att
          mh.z      <- mh.update$mh

          mh.update <- mhupdate(acc=acc.phi.z, att=att.phi.z, mh=mh.phi.z)
          acc.phi.z <- mh.update$acc
          att.phi.z <- mh.update$att
          mh.phi.z  <- mh.update$mh
        }

      }  # fi temporalz
    }  # fi skew

  }  # end nthin

  if (iter > burn) {
    mu <- x.beta + lambda.1* zg
    res <- y - mu
    # predictions
    if (predictions) {
    	yp <- predictY(d11=d11, d12=d12, cov.model=cov.model, rho=rho, nu=nu,
                     gamma=gamma, res=res, beta=beta, tau=tau, taug=taug, z=z,
                     prec=prec, lambda.1=lambda.1, s.pred=s.pred, x.pred=x.pred,
                     knots=knots)

    }
  }

  # storage
  keepers.tau[iter, , ]   <- tau
  keepers.beta[iter, ]    <- beta
  keepers.tau.alpha[iter] <- tau.alpha
  keepers.tau.beta[iter]  <- tau.beta
  keepers.rho[iter]       <- rho
  keepers.nu[iter]        <- nu
  keepers.gamma[iter]     <- gamma
  if (keep.knots & (nknots > 1)) {
    keepers.knots[iter, , , ] <- knots
  }

  if (predictions & iter > burn) {
    y.pred[(iter - burn), , ] <- yp
  }
  if (skew) {
    keepers.lambda.1[iter] <- lambda.1
    keepers.lambda.2[iter] <- lambda.2
    keepers.lambda[iter] <- lambda
    keepers.z[iter, , ]   <- z
  }
  if ((nknots > 1) & !fixknots) {
    keepers.avgparts[iter, ] <- avgparts
  }
  if (temporalw) {
  	keepers.phi.w[iter] <- phi.w
  }
  if (temporalz) {
  	keepers.phi.z[iter] <- phi.z
  }
  if (temporaltau) {
  	keepers.phi.tau[iter] <- phi.tau
  }

  # update notifications for printing
  if (iter %% update == 0) {
  	if (temporalw) {
  	  acc.rate.phi.w <- round(acc.phi.w / att.phi.w, 3)
  	}
  	if (temporalz) {
  	  acc.rate.phi.z <- round(acc.phi.z / att.phi.z, 3)
  	  acc.rate.z <- round(acc.z / att.z, 3)
  	}
   	if (temporaltau) {
   	  acc.rate.phi.tau <- round(acc.phi.tau / att.phi.tau, 3)
      acc.rate.tau     <- round(acc.tau / att.tau, 3)
   	} else {
      acc.rate.tau <- round(acc.tau / att.tau, 3)
    }

  	acc.rate.rho <- round(acc.rho / att.rho, 3)
  	acc.rate.nu <- round(acc.nu / att.nu, 3)
  	acc.rate.gamma <- round(acc.gamma / att.gamma, 3)

  	if (iter < burn) {
  	  begin <- max(1, (iter - 2000))
  	} else {
  	  begin <- burn
  	}

    # what all should be plotted (at most 18 plots)
    # always: beta, rho, nu, gamma, tau.alpha, tau.beta and a few tau terms
    # if multiple knots: 2 knots worth
    # if ts: phi.w, phi.tau, phi.z
    # if skew: lambda, a few z terms
    if (iterplot) {
      if (skew) {
        par(mfrow=c(3, 6))
      } else {
        par(mfrow=c(2, 6))
      }

      plot(keepers.beta[begin:iter, 1], type="l", main="beta 0",
           xlab="", ylab="")
      if (p > 1) {
        plot(keepers.beta[begin:iter, 2], type="l", main="beta 1",
           xlab="", ylab="")
      }
      if (p > 2) {
        plot(keepers.beta[begin:iter, 3], type="l", main="beta 2",
           xlab="", ylab="")
      }

      ylab.rho <- paste("acc =", acc.rate.rho)
      if (mh.rho < 0.00001) {
        xlab.rho <- "mh < 0.00001"
      } else if (mh.rho < 10000) {
        xlab.rho <- paste("mh =", round(mh.rho, 5))
      } else {
        xlab.rho <- "mn > 10000"
      }
      plot(keepers.rho[begin:iter], type="l", main="rho",
           xlab=xlab.rho, ylab=ylab.rho)

      ylab.nu <- paste("acc =", acc.rate.nu)
      if (mh.nu < 0.00001) {
        xlab.nu <- "mh < 0.00001"
      } else if (mh.nu < 10000) {
        xlab.nu <- paste("mh =", round(mh.nu, 5))
      } else {
        xlab.nu <- "mh > 10000"
      }
      plot(keepers.nu[begin:iter], type="l", main="nu",
           xlab=xlab.nu, ylab=ylab.nu)

      ylab.gamma <- paste("acc =", acc.rate.gamma)
      if (mh.gamma < 0.00001) {
        xlab.gamma <- "mh < 0.0001"
      } else if (mh.gamma < 10000) {
        xlab.gamma <- paste("mh =", round(mh.gamma, 5))
      } else {
        xlab.gamma <- "mh > 10000"
      }
      plot(keepers.gamma[begin:iter], type="l", main="gamma",
           xlab=xlab.gamma, ylab=ylab.gamma)

      nparts.1 <- length(which(g[, 1] == 1))
      nparts.2 <- length(which(g[, 10] == 1))
      nparts.3 <- length(which(g[, 21] == 1))
      mh.disp.1 <- round(mh.tau[(nparts.1 + 1)], 2)
      mh.disp.2 <- round(mh.tau[(nparts.2 + 1)], 2)
      mh.disp.3 <- round(mh.tau[(nparts.3 + 1)], 2)
      ylab.tau.1 <- paste("acc = ", acc.rate.tau[1, 1])
      ylab.tau.2 <- paste("acc = ", acc.rate.tau[1, 10])
      ylab.tau.3 <- paste("acc = ", acc.rate.tau[1, 21])

      plot(tau[1, ], type="l", main="tau 1 (all days)")
      if (nknots >= 3) {
        plot(tau[2, ], type="l", main="tau 2 (all days)")
        plot(tau[3, ], type="l", main="tau 3 (all days)")
      }
      # plot(keepers.tau[begin:iter, 1, 1], type="l", main="tau 1,1",
      #      xlab=paste(nparts.1, ", ", mh.disp.1), ylab=ylab.tau.1)
      # plot(keepers.tau[begin:iter, 1, 10], type="l", main="tau 1, 10",
      #      xlab=paste(nparts.2, ", ", mh.disp.2), ylab=ylab.tau.2)
      # plot(keepers.tau[begin:iter, 1, 21], type="l", main="tau 1, 21",
      #      xlab=paste(nparts.3, ", ", mh.disp.3), ylab=ylab.tau.3)

      plot(keepers.tau.alpha[begin:iter], type="l", main="tau.alpha",
           xlab="", ylab="")
      plot(keepers.tau.beta[begin:iter], type="l", main="tau.beta",
           xlab="", ylab="")

      if (skew) {
        plot(keepers.lambda[begin:iter], type="l", main="lambda",
             xlab="", ylab="")

        if (temporalz) {
          ylab.z.1 <- paste("acc =", acc.rate.z[1, 1])
          ylab.z.2 <- paste("acc =", acc.rate.z[1, 10])
          ylab.z.3 <- paste("acc =", acc.rate.z[1, 21])
        } else {
          ylab.z.1 <- ""
          ylab.z.2 <- ""
          ylab.z.3 <- ""
        }
        plot(z[1, ], type="l", main="z 1 (all days)")
        if (nknots >= 3) {
          plot(z[2, ], type="l", main="z 2 (all days)")
          plot(z[3, ], type="l", main="z 3 (all days)")
        }
        # plot(keepers.z[begin:iter, 1, 1], type="l", main="z 1, 1",
        #      xlab="", ylab=ylab.z.1)
        # plot(keepers.z[begin:iter, 1, 10], type="l", main="z 1, 10",
        #      xlab="", ylab=ylab.z.2)
        # plot(keepers.z[begin:iter, 1, 21], type="l", main="z 1, 21",
        #      xlab="", ylab=ylab.z.3)
      }

      # if (temporalw) {
      #   ylab.phi.w <- paste("acc =", acc.rate.phi.w)
      #   plot(keepers.phi.w[begin:iter], type="l",
      #        main="phi.w",
      #        xlab=paste("mh =", round(mh.phi.w, 3)), ylab=ylab.phi.w)
      # }
      # if (temporalz) {
      #   ylab.phi.z <- paste("acc =", acc.rate.phi.z)
      # 	plot(keepers.phi.z[begin:iter], type="l",
      #        main="phi.z",
      #        xlab=paste("mh =", round(mh.phi.z, 3)), ylab=ylab.phi.z)
      # }
      # if (temporaltau) {
      #   ylab.phi.tau <- paste("acc =", acc.rate.phi.tau)
      # 	plot(keepers.phi.tau[begin:iter], type="l",
      #        main="phi.tau",
      #        xlab=paste("mh =", round(mh.phi.tau, 3)), ylab=ylab.phi.tau)
      # }
      hist(y.init, main="true y")
      hist(y, main="imputed y")

    }

    cat("\t iter", iter, "\n")
  }


  } #end iters

if ((nknots == 1) | fixknots) {
  keepers.avgparts <- NULL
} else {
  keepers.avgparts <- keepers.avgparts[return.iters, ]
}

if (keep.knots & (nknots > 1)) {
  keepers.knots <- keepers.knots[return.iters, , , ]
} else {
  keepers.knots <- NULL
}

if (!predictions) {
  y.pred <- NULL
}

if (!skew) {
  keepers.z <- NULL
  keepers.lambda <- NULL
} else {
  keepers.z <- keepers.z[return.iters, , ]
  keepers.lambda <- keepers.lambda[return.iters]
}

if (!temporalz) {  # ts
  keepers.phi.z <- NULL
} else {
  keepers.phi.z <- keepers.phi.z[return.iters]
}
if (!temporalw) {  # ts
  keepers.phi.w <- NULL
} else {
  keepers.phi.w <- keepers.phi.w[return.iters]
}
if (!temporaltau) {  # ts
  keepers.phi.tau <- NULL
} else {
  keepers.phi.tau <- keepers.phi.tau[return.iters]
}

results <- list(tau=keepers.tau[return.iters, , ],
                beta=keepers.beta[return.iters, ],
                tau.alpha=keepers.tau.alpha[return.iters],
                tau.beta=keepers.tau.beta[return.iters],
                rho=keepers.rho[return.iters],
                nu=keepers.nu[return.iters],
                gamma=keepers.gamma[return.iters],
                yp=y.pred,
                lambda=keepers.lambda,
                z=keepers.z,
                knots=keepers.knots,
                avgparts=keepers.avgparts,
                phi.z=keepers.phi.z,  # ts
                phi.w=keepers.phi.w,  # ts
                phi.tau=keepers.phi.tau  # ts
                )

return(results)
}#end mcmc()

