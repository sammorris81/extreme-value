################################################################################
# MCMC
#
# Let ns be the number of sites, nt be the number of days, np be the number
# of parameters, and npred be the number of sites at which to make predictions.
# Finally, the model is written using sigma^2 as variance, but from a
# computation perspective, it's easier to store and use tau.
#
# y(ns, nt): data
# s(ns, 2): sites
# x(ns, nt, np): covariates
# s.pred(npred, 2): sites at which to make predictions
# x.pred(npred, nt, np): covariates for sites at which to make predictions
# min.s(2): smallest x and y
# max.s(2): largest x and y
# thresh.all(1): a scalar for value of the threshold (thresh.quant = FALSE) or
#                the quantile for the threshold (thresh.quant = TRUE)
# thresh.quant(1): boolean for whether thresholds given are quantiles (0, 1) or
#                  actual values from the data
# thresh.site.specific(1): boolean for site-specific threshold
# thresh.site(ns): vector of site-specific thresholds
# nknots(1): the number of knots to use in the MCMC
# fixknots(1): boolean on whether knot locations should be fixed
# keep.knots(1): should the MCMC return the knots in the output
# skew(1): boolean for whether a skew-model should be fit
# method(1): string ("t" or "gaussian") to indicate which process to fit
# cov.model(1): string ("matern" or "exponential")
# temporalw(1): boolean for TS on the knot location
# temporaltau(1): boolean for TS on the precision
# temporalz(1): boolean for TS on the skew term
#
# Starting values:
# beta.init(np): vector of starting values for beta
# tau.init(1 or nt): scalar or vector of starting values for precision
# tau.alpha.init(1): starting value of the a hyperparameter for precision
# tau.beta.init(1): starting value of the b hyperparameter for precision
# rho.init(1): starting value of the bandwidth for Matern Sigma
# nu.init(1): starting value of the smoothness for Matern Sigma
# gamma.init(1): starting value for strength of the nugget effect
# z.init(1 or nknots): starting values for z
# lambda.init(1): starting value for lambda
#
# Prior distributions:
# beta      ~iid N(beta.m, beta.s)
# tau.alpha ~ discrete(seq(tau.alpha.min, tau.alpha.max, by = tau.alpha.by))
# tau.beta  ~ gamma(tau.beta.a, tau.beta.b)
# rho       ~ uniform(0, rho.upper)
# lambda    ~ N(lambda.m, lambda.s)
#
# MCMC Settings:
# iters(1): number of iterations to run the MCMC
# burn(1): length of burnin period
# update(1): how often to update MCMC progress
# thin(1): how many observations to thin in the MCMC
# iterplot(1): boolean for whether iteration plots should be generated at an
#              update iteration for the MCMC
#
################################################################################
if (!exists("update_params_cpp_init")) {
  source('update_params_cpp.R')
}
source('update_params.R')

mcmc <- function(y, s, x, s.pred = NULL, x.pred = NULL,
                 min.s, max.s,  # don't want to specify defaults
                 thresh.all = 0, thresh.quant = TRUE,
                 thresh.site.specific = FALSE, thresh.site = NULL,
                 nknots = 1, fixknots = FALSE,
                 keep.knots = FALSE,
                 skew = T, method = "t",  # or "gaussian"
                 cov.model = "matern",  # or "exponential"
                 temporalw = TRUE, temporaltau = TRUE, temporalz = TRUE,
                 # initial values
                 beta.init = NULL, tau.init = 1,
                 tau.alpha.init = 0.1, tau.beta.init = 0.1,
                 rho.init = 1, nu.init = 0.5, gamma.init = 0.5,
                 z.init = 0, lambda.init = NULL, knots.init = NULL,
                 # priors
                 beta.m = 0, beta.s = 20,
                 tau.alpha.min = 0.2, tau.alpha.max = 20, tau.alpha.by = 0.2,
                 tau.beta.a = 1, tau.beta.b = 1,
                 rho.upper = NULL,
                 lognu.m = -1.2, lognu.s = 1, nu.upper = NULL, fix.nu = FALSE,
                 lambda.m = 0, lambda.s = 20,
                 # mcmc settings
                 iters = 5000, burn = 1000, update = 100, thin = 1,
                 iterplot = FALSE
){

  library(SpatialTools)
  library(fields)
  library(emulator)
  library(e1071)  # for skewness function to get initial lambda if not set
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
    y.pred <- array(0, c((iters - burn), np, nt))
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
      warning("Warning: setting site-specific thresholds to thresh.all")
      thresh.site <- thresh.all
    }
    if (thresh.quant) {
      if ((thresh.site < 0) | (thresh.site > 1)) {
        stop("Error: thresh.site should be the desired site-specific quantile
            between 0 and 1")
      }
      # if it's site specific, then we want to keep everything over the 95th
      # quantile for the data set and the max for each site in a matrix that's
      # ns x nt
      thresh.site <- apply(y, 1, quantile, probs = thresh.site, na.rm = T)
    }
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

  missing.obs <- is.na(y)
  if (sum(missing.obs) > 0) {
    y[missing.obs]  <- mean(y, na.rm = TRUE)
    missing <- TRUE
  } else {
    missing <- FALSE
  }

  # initialize partition
  if (!fixknots) {
    if (is.null(knots.init)) {
      # constrain the initial draw to be near the middle for stability
      lower.1 <- min.s[1] + (max.s[1] - min.s[1]) / 6
      upper.1 <- max.s[1] - (max.s[1] - min.s[1]) / 6
      lower.2 <- min.s[2] + (max.s[2] - min.s[2]) / 6
      upper.2 <- max.s[2] - (max.s[2] - min.s[2]) / 6
      knots <- array(NA, dim = c(nknots, 2, nt))
      knots[, 1, ] <- runif(nknots * nt, lower.1, upper.1)
      knots[, 2, ] <- runif(nknots * nt, lower.2, upper.2)
    } else {
      knots <- array(knots.init, dim = c(nknots, 2, nt))
    }
    knots.star        <- array(NA, c(nknots, 2, nt))
    knots.star[, 1, ] <- transform$probit(knots[, 1, ], lower = min.s[1],
                                          upper = max.s[1])
    knots.star[, 2, ] <- transform$probit(knots[, 2, ], lower = min.s[2],
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
    beta[1] <- mean(y)
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

  zg <- matrix(0, ns, nt)
  if (length(z.init) == 1 && skew) {
    cat("\t initializing all z terms to", z.init, "\n")
  }
  z <- matrix(z.init, nknots, nt)
  for (t in 1:nt) {
    zg[, t] <- z[g[, t], t]
  }

  if (skew) {
    if (is.null(lambda.init)) {
      lambda <- skewness(y, na.rm = TRUE)
    } else {
      lambda <- lambda.init
    }
  } else {
    lambda <- 0
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

  fixnu <- FALSE
  if (cov.model == "exponential") {
    nu <- 0.5
    fix.nu <- TRUE
  }

  #### Precomputed initialized values.
  # setup correlation matrix without multiplying by gamma
  # simple.cov.sp only needs to be updated for rho and nu updates
  cor <- simple.cov.sp(D = d, sp.type = "matern", sp.par = c(1, rho),
                       error.var = 0, smoothness = nu, finescale.var = 0)
  C <- CorFx(d = d, gamma = gamma, rho = rho, nu = nu)
  CC <- tryCatch(chol.inv(C, inv = TRUE, logdet = TRUE),
                 error = function(e) {
                   eig.inv(C, inv = TRUE, logdet = TRUE, mtx.sqrt = TRUE)
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

  acc.tau <- att.tau <- matrix(1, nrow = nknots, ncol = nt)
  mh.tau <- matrix(0.05, nknots, nt)
  nparts.tau <- matrix(1, nrow = nknots, ncol = nt)

  if (temporalz) {
    phi.z <- 0
    acc.z <- att.z <- mh.z <- matrix(1, nknots, nt)
    acc.phi.z <- att.phi.z <- mh.phi.z <- 1
  }

  # MH tuning params
  acc.w      <- att.w      <- mh.w     <- matrix(0.1, nknots, nt)  # knot locs
  acc.delta  <- att.delta  <- mh.delta <- 0.1
  acc.rho    <- att.rho    <- mh.rho   <- 1
  acc.nu     <- att.nu     <- mh.nu    <- 0.1
  acc.gamma  <- att.gamma  <- mh.gamma <- 1

  # storage
  keepers.tau       <- array(NA, dim = c(iters, nknots, nt))
  keepers.beta      <- matrix(NA, nrow = iters, ncol = p)
  keepers.ll        <- matrix(NA, nrow = iters, ncol = nt)
  keepers.tau.alpha <- rep(NA, iters)
  keepers.tau.beta  <- rep(NA, iters)
  keepers.nparts    <- matrix(NA, nrow = iters, ncol = nt)
  keepers.rho       <- rep(NA, iters)
  keepers.nu        <- rep(NA, iters)
  keepers.gamma     <- rep(NA, iters)
  keepers.z         <- array(NA, dim=c(iters, nknots, nt))
  keepers.lambda    <- rep(NA, iters)

  keepers.avgparts  <- matrix(NA, nrow = iters, ncol = nt)  # avg parts per day
  if (keep.knots & (nknots > 1)) {
    keepers.knots     <- array(NA, dim = c(iters, nknots, 2, nt))
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
  keepers.y    <- array(0, c((iters - burn), ns, nt))
  return.iters <- (burn + 1):iters

  tic <- proc.time()
  for (iter in 1:iters) { for (ttt in 1:thin) {

    # data imputation
    if (thresholded) {  # do data imputation and store as y
      mu <- x.beta + lambda * zg
      res <- y - mu
      y <- imputeY(y = y, taug = taug, mu = mu, obs = thresh.obs, cor = cor,
                   gamma = gamma, thresh.mtx = thresh.mtx)
    }

    # missing values
    if (missing) {
      mu <- x.beta + lambda * zg
      res <- y - mu
      y <- imputeY(y = y, taug = taug, mu = mu, obs = missing.obs, cor = cor,
                   gamma = gamma)
    }

    # update beta
    mu <- x.beta + lambda * zg
    res <- y - mu
    this.update <- updateBeta(beta.m = beta.m, beta.s = beta.s,
                              x = x, y = y, zg = zg, taug = taug,
                              prec = prec, skew = skew,
                              lambda = lambda, lambda.m = lambda.m,
                              lambda.s = lambda.s)
    beta   <- this.update$beta
    lambda <- this.update$lambda

    for (t in 1:nt) {
      x.beta[, t] <- x[, t, ] %*% beta
    }
    mu <- x.beta + lambda * zg
    res <- y - mu

    # update partitions
    if ((nknots > 1) & (!fixknots)) {
      avgparts <- rep(0, nt)
      this.update <- updateKnotsTS(phi = phi.w, knots = knots,
                                   g = g, ts = temporalw,
                                   tau = tau, z = z, s = s,
                                   min.s = min.s, max.s = max.s,
                                   x.beta = x.beta,
                                   lambda = lambda, y = y,
                                   prec = prec, att = att.w,
                                   acc = acc.w, mh = mh.w,
                                   att.phi = att.phi.w,
                                   acc.phi = acc.phi.w,
                                   mh.phi = mh.phi.w)
      knots.star <- this.update$knots.star
      knots      <- this.update$knots
      g          <- this.update$g
      taug       <- this.update$taug
      zg         <- this.update$zg
      acc.w      <- this.update$acc
      att.w      <- this.update$att
      phi.w      <- this.update$phi
      acc.phi.w  <- this.update$acc.phi
      att.phi.w  <- this.update$att.phi

      # update metropolis sds
      if (iter < (burn / 2)) {
        this.update  <- mhupdate(acc = acc.w, att = att.w, mh = mh.w,
                                 nattempts = nknots * nt)
        acc.w        <- this.update$acc
        att.w        <- this.update$att
        mh.w         <- this.update$mh

        if (temporalw) {
          this.update <- mhupdate(acc = acc.phi.w, att = att.phi.w, mh = mh.phi.w)
          acc.phi.w   <- this.update$acc
          att.phi.w   <- this.update$att
          mh.phi.w    <- this.update$mh
        }
      }
    }  # fi nknots > 1

    # update tau
    # all taus require mu and res
    mu <- x.beta + lambda * zg
    res <- y - mu
    if (method == "gaussian") {  # single random effect for all days
      this.update <- updateTauGaus(res = res, prec = prec, tau.alpha = tau.alpha,
                                   tau.beta = tau.beta)
      tau  <- matrix(this.update, nknots, nt)
      taug <- matrix(this.update, ns, nt)

    } else if (method == "t") {
      if (!temporaltau) {
        this.update <- updateTau(tau = tau, taug = taug, g = g,
                                 res = res, nparts.tau = nparts.tau,
                                 prec = prec, z = z,
                                 tau.alpha = tau.alpha,
                                 tau.beta = tau.beta, skew = skew,
                                 att = att.tau, acc = acc.tau,
                                 mh = mh.tau)
        tau     <- this.update$tau
        taug    <- this.update$taug
        acc.tau <- this.update$acc
        att.tau <- this.update$att

        # # Note: this doesn't really tune the algorithm because mh tau isn't used
        # if (iter < (burn / 2)) {
        #   mh.update  <- mhupdate(acc=acc.tau, att=att.tau, mh=mh.tau)
        #   acc.tau <- mh.update$acc
        #   att.tau <- mh.update$att
        #   mh.tau  <- mh.update$mh
        # }
      } else {
        this.update <- updateTauTS(phi = phi.tau, tau = tau,
                                   taug = taug, g = g, res = res,
                                   nparts.tau = nparts.tau,
                                   prec = prec, z = z,
                                   tau.alpha = tau.alpha,
                                   tau.beta = tau.beta, skew = skew,
                                   att = att.tau, acc = acc.tau,
                                   mh = mh.tau,
                                   att.phi = att.phi.tau,
                                   acc.phi = acc.phi.tau,
                                   mh.phi = mh.phi.tau)

        tau         <- this.update$tau
        taug        <- this.update$taug
        phi.tau     <- this.update$phi
        acc.tau     <- this.update$acc
        att.tau     <- this.update$att
        acc.phi.tau <- this.update$acc.phi
        att.phi.tau <- this.update$att.phi


        # update metropolis sds
        if (iter < (burn / 2)) {
          this.update <- mhupdate(acc = acc.tau, att = att.tau, mh = mh.tau)
          acc.tau     <- this.update$acc
          att.tau     <- this.update$att
          mh.tau      <- this.update$mh

          this.update  <- mhupdate(acc = acc.phi.tau, att = att.phi.tau,
                                   mh = mh.phi.tau)
          acc.phi.tau <- this.update$acc
          att.phi.tau <- this.update$att
          mh.phi.tau  <- this.update$mh
        }
      }

      # update tau.alpha and tau.beta
      tau.alpha <- updateTauAlpha(tau = tau, tau.beta = tau.beta,
                                  tau.alpha.min = tau.alpha.min,
                                  tau.alpha.max = tau.alpha.max,
                                  tau.alpha.by = tau.alpha.by)
      tau.beta  <- updateTauBeta(tau = tau, tau.alpha = tau.alpha,
                                 tau.beta.a = tau.beta.a,
                                 tau.beta.b = tau.beta.b)
    }  # fi method == t

    # update rho and nu and gamma
    # update rss with new residuals and taug
    mu <- x.beta + lambda * zg
    res <- y - mu
    cur.rss <- sum(rss(prec = prec, y = sqrt(taug) * res))

    # update rho
    this.update <- updateRho(rho = rho, nu = nu, d = d, gamma = gamma,
                             res = res, taug = taug, prec = prec, cor = cor,
                             logdet.prec = logdet.prec, cur.rss = cur.rss,
                             rho.upper = rho.upper,
                             att = att.rho, acc = acc.rho, mh = mh.rho)
    rho         <- this.update$rho
    prec        <- this.update$prec
    cor         <- this.update$cor
    logdet.prec <- this.update$logdet.prec
    cur.rss     <- this.update$cur.rss
    acc.rho     <- this.update$acc
    att.rho     <- this.update$att

    # update nu
    if (!fix.nu) {
      this.update <- updateNu(nu = nu, lognu.m = lognu.m, lognu.s = lognu.s,
                              rho = rho, d = d, gamma = gamma, res = res,
                              taug = taug, prec = prec, cor = cor,
                              logdet.prec = logdet.prec, cur.rss = cur.rss,
                              nu.upper = nu.upper,
                              att = att.nu, acc = acc.nu, mh = mh.nu)
      nu          <- this.update$nu
      prec        <- this.update$prec
      cor         <- this.update$cor
      logdet.prec <- this.update$logdet.prec
      cur.rss     <- this.update$cur.rss
      att.nu      <- this.update$att
      acc.nu      <- this.update$acc
    }


    # gamma - make sure the cor argument does not include gamma
    this.update <- updateGamma(gamma = gamma, d = d, rho = rho, nu = nu,
                               taug = taug, res = res, prec = prec, cor = cor,
                               logdet.prec = logdet.prec, cur.rss = cur.rss,
                               att = att.gamma, acc = acc.gamma, mh = mh.gamma)
    gamma       <- this.update$gamma
    prec        <- this.update$prec
    logdet.prec <- this.update$logdet.prec
    cur.rss     <- this.update$cur.rss
    acc.gamma   <- this.update$acc
    att.gamma   <- this.update$att

    if (iter < (burn / 2)) {
      this.update <- mhupdate(acc = acc.rho, att = att.rho, mh = mh.rho)
      acc.rho     <- this.update$acc
      att.rho     <- this.update$att
      mh.rho      <- this.update$mh

      this.update <- mhupdate(acc = acc.nu, att = att.nu, mh = mh.nu)
      acc.nu      <- this.update$acc
      att.nu      <- this.update$att
      mh.nu       <- this.update$mh

      this.update <- mhupdate(acc = acc.gamma, att = att.gamma, mh = mh.gamma)
      acc.gamma   <- this.update$acc
      att.gamma   <- this.update$att
      mh.gamma    <- this.update$mh
    }

    # update z (skew random effects)
    # Note: lambda is updated with beta terms
    if (skew) {
      mu <- x.beta + lambda * zg
      res <- y - mu

      if (!temporalz) {
        this.update <- updateZ(y = y, x.beta = x.beta, zg = zg,
                               prec = prec, tau = tau, mu = mu,
                               taug = taug, g = g, lambda = lambda)

        z  <- this.update$z
        zg <- this.update$zg
      } else {
        # trying to do the update with lambda = lambda.1 * lambda.2
        this.update <- updateZTS(z = z, zg = zg, y = y,
                                 lambda = lambda, x.beta = x.beta,
                                 phi = phi.z, tau = tau,
                                 taug = taug, g = g, prec = prec,
                                 acc = acc.z, att = att.z,
                                 mh = mh.z, acc.phi = acc.phi.z,
                                 att.phi = att.phi.z,
                                 mh.phi = mh.phi.z)
        z      <- this.update$z
        zg     <- this.update$zg
        att.z  <- this.update$att
        acc.z  <- this.update$acc
        phi.z  <- this.update$phi
        att.phi.z <- this.update$att.phi
        acc.phi.z <- this.update$acc.phi

        if (iter < (burn / 2)) {
          this.update <- mhupdate(acc = acc.z, att = att.z, mh = mh.z)
          acc.z       <- this.update$acc
          att.z       <- this.update$att
          mh.z        <- this.update$mh

          this.update <- mhupdate(acc = acc.phi.z, att = att.phi.z,
                                  mh = mh.phi.z)
          acc.phi.z   <- this.update$acc
          att.phi.z   <- this.update$att
          mh.phi.z    <- this.update$mh
        }

      }  # fi temporalz
    }  # fi skew

  }  # end nthin

    if (iter > burn) {
      # predictions
      if (predictions) {
        mu <- x.beta + lambda * zg
        res <- y - mu
        yp <- predictY(d11 = d11, d12 = d12, cov.model = cov.model,
                       rho = rho, nu = nu, gamma = gamma, res = res,
                       beta = beta, tau = tau, taug = taug, z = z,
                       prec = prec, lambda = lambda, s.pred = s.pred,
                       x.pred = x.pred, knots = knots)

      }
    }

    # storage
    keepers.tau[iter, , ]   <- tau
    keepers.beta[iter, ]    <- beta
    keepers.tau.alpha[iter] <- tau.alpha * 2  # in the manuscript use a / 2
    keepers.tau.beta[iter]  <- tau.beta * 2   # in the manuscript use b / 2
    keepers.rho[iter]       <- rho
    keepers.nu[iter]        <- nu
    keepers.gamma[iter]     <- gamma
    if (keep.knots & (nknots > 1)) {
      keepers.knots[iter, , , ] <- knots
    }

    if (predictions & iter > burn) {
      y.pred[(iter - burn), , ] <- yp
    }
    if (iter > burn) {
      keepers.y[(iter - burn), , ] <- y
    }
    if (skew) {
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
          par(mfrow = c(3, 6))
        } else {
          par(mfrow = c(2, 6))
        }

        plot(keepers.beta[begin:iter, 1], type = "l", main = "beta 0",
             xlab = "", ylab = "")
        if (p > 1) {
          plot(keepers.beta[begin:iter, 2], type = "l", main = "beta 1",
               xlab = "", ylab = "")
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
        plot(keepers.rho[begin:iter], type = "l", main = "rho",
             xlab = xlab.rho, ylab = ylab.rho)

        ylab.nu <- paste("acc =", acc.rate.nu)
        if (mh.nu < 0.00001) {
          xlab.nu <- "mh < 0.00001"
        } else if (mh.nu < 10000) {
          xlab.nu <- paste("mh =", round(mh.nu, 5))
        } else {
          xlab.nu <- "mh > 10000"
        }
        plot(keepers.nu[begin:iter], type = "l", main = "nu",
             xlab = xlab.nu, ylab = ylab.nu)

        ylab.gamma <- paste("acc =", acc.rate.gamma)
        if (mh.gamma < 0.00001) {
          xlab.gamma <- "mh < 0.0001"
        } else if (mh.gamma < 10000) {
          xlab.gamma <- paste("mh =", round(mh.gamma, 5))
        } else {
          xlab.gamma <- "mh > 10000"
        }
        plot(keepers.gamma[begin:iter], type = "l", main = "gamma",
             xlab = xlab.gamma, ylab = ylab.gamma)

        if (nt <= 10) {
          nparts.1 <- length(which(g[, 1] == 1))
          mh.disp.1 <- round(mh.tau[(nparts.1 + 1)], 2)
          ylab.tau.1 <- paste("acc = ", acc.rate.tau[1, 1])
        } else if (nt <= 20) {
          nparts.2 <- length(which(g[, 11] == 1))
          mh.disp.2 <- round(mh.tau[(nparts.2 + 1)], 2)
          ylab.tau.2 <- paste("acc = ", acc.rate.tau[1, 11])
        } else {
          nparts.3 <- length(which(g[, 21] == 1))
          mh.disp.3 <- round(mh.tau[(nparts.3 + 1)], 2)
          ylab.tau.3 <- paste("acc = ", acc.rate.tau[1, 21])
        }

        plot(tau[1, ], type="l", main = "tau 1 (all days)")
        if (nknots >= 3) {
          plot(tau[2, ], type="l", main = "tau 2 (all days)")
        }

        plot(keepers.tau.alpha[begin:iter], type = "l", main = "tau.alpha",
             xlab = "", ylab = "")
        plot(keepers.tau.beta[begin:iter], type = "l", main = "tau.beta",
             xlab = "", ylab = "")

        if (skew) {
          plot(keepers.lambda[begin:iter], type = "l", main = "lambda",
               xlab = "", ylab = "")

          if (temporalz) {
            if (nt <= 10) {
              ylab.z.1 <- paste("acc =", acc.rate.z[1, 1])
            } else if (nt <= 20) {
              ylab.z.2 <- paste("acc =", acc.rate.z[1, 11])
            } else {
              ylab.z.3 <- paste("acc =", acc.rate.z[1, 21])
            }
          } else {
            ylab.z.1 <- ""
            ylab.z.2 <- ""
            ylab.z.3 <- ""
          }
          plot(z[1, ], type = "l", main = "z 1 (all days)")
          if (nknots >= 3) {
            plot(z[2, ], type = "l", main = "z 2 (all days)")
          }
        }

        if (temporalz) {
          ylab.phi.z <- paste("acc =", acc.rate.phi.z)
          plot(keepers.phi.z[begin:iter], type = "l", main = "phi.z",
               xlab = paste("mh =", round(mh.phi.z, 3)), ylab = ylab.phi.z)
        }
        if (temporaltau) {
          ylab.phi.tau <- paste("acc =", acc.rate.phi.tau)
          plot(keepers.phi.tau[begin:iter], type = "l", main = "phi.tau",
               xlab=paste("mh =", round(mh.phi.tau, 3)), ylab = ylab.phi.tau)
        }
        hist(y.init, main = "true y")
        hist(y, main = "imputed y")

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

  results <- list(tau = keepers.tau[return.iters, , ],
                  beta = keepers.beta[return.iters, ],
                  tau.alpha = keepers.tau.alpha[return.iters],
                  tau.beta = keepers.tau.beta[return.iters],
                  rho = keepers.rho[return.iters],
                  nu = keepers.nu[return.iters],
                  gamma = keepers.gamma[return.iters],
                  y = keepers.y,
                  yp = y.pred,
                  lambda = keepers.lambda,
                  z = keepers.z,
                  knots = keepers.knots,
                  avgparts = keepers.avgparts,
                  phi.z = keepers.phi.z,  # ts
                  phi.w = keepers.phi.w,  # ts
                  phi.tau = keepers.phi.tau  # ts
  )

  return(results)
}  # end mcmc()

