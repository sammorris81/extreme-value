#########################################################################
# MCMC 
#
# TODO: Add in model description here
#
#ASSUMES x and y are in [0,1]^2
#########################################################################
source('condmean_cpp.R')

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL, 
                 thresh.all=0, thresh.quant=T, nknots=1, keep.knots=F,
                 iters=5000, burn=1000, update=100, thin=1,
                 iterplot=F, plotname=NULL, method="t", 
                 # just to debug temporal parts.
                 temporalw=F, temporaltau=F, temporalz=F,  # eventually change to temporal=F
                 # initial values
                 beta.init=NULL, tau.init=2, tau.alpha.init=0.1, tau.beta.init=0.1,
                 rho.init=5, nu.init=0.5, gamma.init=0.5,
                 # priors
                 beta.m=0, beta.s=10, 
                 tau.alpha.m=0, tau.alpha.s=1, 
                 tau.beta.a=0.1, tau.beta.b=0.1,
                 logrho.m=0, logrho.s=10, rho.upper=NULL,
                 lognu.m=-1.2, lognu.s=1, nu.upper=NULL,
                 gamma.m=0, gamma.s=1,
                 # covariance model
                 cov.model="matern",  # or "exponential"
                 rho.prior="cont",  # or "disc"
                 # skew inits
                 z.init=1, lambda.init=0,
                 # skew priors
                 lambda.m=0, lambda.s=10, skew=T,
                 thresh.site.specific=F, thresh.site=NULL,
                 # troubleshooting
                 debug=F, fixhyper=F, tau.t, z.t
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
    y.pred <- array(0, c(iters, np, nt))
  }

  d       <- rdist(s)  # distance between sites
  diag(d) <- 0
    
  # store the loc/day for observed values below thresh.
  if (thresh.quant & thresh.all > 0) {  # threshold based on sample quantiles
    thresh.all.q <- quantile(y, thresh.all, na.rm=T)
  } else {
    thresh.all.q <- thresh.all
  }
  thresh.all.mtx <- matrix(thresh.all.q, ns, nt)  # want as a matrix for easy replacement
  
  if (thresh.site.specific) { 
  	if (is.null(thresh.site)) {
  	  warning("Warning: setting site-specific time series threshold to thresh.")
  	  thresh.site <- thresh.all
  	}
  	if ((thresh.site < 0) | (thresh.site > 1)) {
  	  stop("Error: thresh.site should be the desired site-specific quantile between 0 and 1")
  	}
    # if it's site specific, then we want to keep everything over the 95th quantile for the data
    # set and the max for each site in a matrix that's ns x nt
    thresh.site <- apply(y, 1, quantile, probs=thresh.site, na.rm=T)
    thresh.site.mtx <- matrix(thresh.site, ns, nt)
    thresh.mtx <- ifelse(
      thresh.site < thresh.all.mtx,
      thresh.site,
      thresh.all.mtx
    )
  } else {  # if the mean doesn't have a site component, use the same threshold for all
  	cat("\t no site-specific threshold set \n")
    thresh.mtx  <- thresh.all.mtx  
  }
  thresh.obs  <- !is.na(y) & (y < thresh.mtx)
  
  missing.obs <- is.na(y)
  y[missing.obs]  <- mean(y, na.rm=T)

    
  # initialize partition
  x.range <- max(s[, 1]) - min(s[, 1])  # gives span of x
  y.range <- max(s[, 2]) - min(s[, 2])  # gives span of y
  range <- max(x.range, y.range)        # want to scale by the same amount in both directions
  knots.con    <- array(rnorm(nknots * nt * 2), c(nknots, 2, nt))
  knots        <- pnorm(knots.con)  # in [0, 1] x [0, 1]
  knots[, 1, ] <- knots[, 1, ] * range + min(s[, 1])  # rescaled back to size of s
  knots[, 2, ] <- knots[, 2, ] * range + min(s[, 2])  # rescaled back to size of s
    
  # initialize parameters
  beta    <- rep(0, p)
  beta[1] <- mean(y)
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
  } else {
  	if (lambda.init != 0) {
  	  warning("lambda.init being ignored since skew=F")
  	}
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
  fixnu <- F
  
  if (cov.model == "exponential") {
    nu <- 0.5
    fixnu <- T
  }
  
  if (rho.prior == "cont") {
    C <- CorFx(d=d, gamma=gamma, rho=rho, nu=nu)
    CC <- tryCatch(chol.inv(C, inv=T, logdet=T),
                   error = function(e) {
                     eig.inv(C, inv=T, logdet=T, mtx.sqrt=T)
                   })
    sd.mtx <- CC$sd.mtx
    prec.cor <- CC$prec
    logdet.prec <- CC$logdet.prec  # already includes 0.5
  } else if (rho.prior == "disc"){  # precompute eigenvectors and eigenvalues
  	if (cov.model != "exponential") {
  	  stop('to use a discrete prior on rho, you must set cov.model = "exponential"')
  	}
  	if (max(s) > 1) {
  		stop('to use a discrete prior on rho, you need to scale your data to a unit square')
  	}
  	
    rhos <- seq(0.01, 1.2, 0.01)  # restricting spatial domain to [0, 1] x [0, 1]
   	C.vectors <- array(NA, dim=c(ns, ns, length(rhos)))
   	C.values  <- matrix(NA, nrow=ns, ncol=length(rhos))
   	for (i in 1:length(rhos)) {
   	  temp.C <- exp(-(d / rhos[i]))
   	  temp.C.eigen <- eigen(temp.C, symmetric=T)
   	  C.vectors[, , i] <- temp.C.eigen$vectors
   	  C.values[, i] <- temp.C.eigen$values
   	}
   	rho.idx <- which(rhos == rho)
   	C <- exp(-(d / rho))
   	D <- 1 - gamma + gamma * C.values[, rho.idx]  # eigenvalues with gamma
   	prec.cor <- quad.tform(diag(1 / D), C.vectors[, , rho.idx])
   	logdet.prec <- -0.5 * sum(log(D))
  } else {
    stop("rho.prior must be cont or disc")
  }
  
  # time series in the random knots/partitions
  if (temporalw) {
    phi.w <- 0
    acc.phi.w <- att.phi.w <- mh.phi.w <- 1
  }
  if (temporaltau) {
    phi.tau <- 0
    acc.phi.tau <- att.phi.tau <- 1
    mh.phi.tau <- 1
    acc.tau.ns <- att.tau.ns <- rep(1, (ns + 1))
    mh.tau.ns <- seq(5, 1, length=(ns + 1))
    acc.tau.low <- acc.tau.high <- 0
  } else {
  	acc.tau.ns <- att.tau.ns <- rep(1, (ns + 1))
  	mh.tau.ns <- seq(5, 1, length=(ns + 1))
    acc.tau.low <- acc.tau.high <- 0
  }
  tau.trials <- nknots * nt * 3
  
  if (temporalz) {
  	z.star <- z  # need a place to keep track of normal values for time series
    phi.z <- 0.5
    acc.z <- att.z <- matrix(1, nknots, nt)
    mh.z <- matrix(15, nknots, nt)
    acc.phi.z <- att.phi.z <- mh.phi.z <- 1
  }
  
  acc.tau <- att.tau   <- matrix(1, nrow=nknots, ncol=nt) 
  mh.tau <- matrix(0.05, nknots, nt)
  nparts.tau <- matrix(1, nrow=nknots, ncol=nt)

  # MH tuning params
  acc.w      <- att.w      <- mh.w     <- rep(0.1, nt)  # knot locations
  acc.delta  <- att.delta  <- mh.delta <- 0.1  
  acc.rho    <- att.rho    <- mh.rho   <- 0.1
  acc.nu     <- att.nu     <- mh.nu    <- 0.1
  acc.gamma  <- att.gamma  <- mh.gamma <- 0.5
  
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
  keepers.lambda   <- rep(NA, iters)
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
    if (thresh.all != 0) {
      mu <- x.beta + lambda * zg
      thresh.mtx.fudge <- 0.99999 * thresh.mtx  # numerical stability
      y.impute <- matrix(y, ns, nt)
      
      for (t in 1:nt) {
      	taug.t <- sqrt(taug[, t])
      	mu.t <- mu[, t]
      	res.t <- y[, t] - mu[, t]
      	impute.these <- which(thresh.obs[, t])
      	
      	# cpp function to find all conditional means and standard deviations
      	impute.cond <- conditional.mean(mn=mu.t, prec=prec.cor, res=res.t, 
      	                                taug=taug.t, include=impute.these)
      	impute.sd <- impute.cond$cond.sd
      	impute.e  <- impute.cond$cond.mn
      	
        u.upper    <- pnorm(thresh.mtx[impute.these, t], impute.e, impute.sd)
        u.impute   <- runif(length(impute.these))
        y.impute.t <- ifelse(  # for numerical stability
          u.upper < 1e-6,
          thresh.mtx.fudge[impute.these, t],
          impute.e + impute.sd * qnorm(u.impute * u.upper)
        )
        y.impute[impute.these, t] <- y.impute.t
        
        # missing values next
        missing.these <- which(missing.obs[, t])  # gives sites that are missing on day t.
        sig.t <- 1 / taug.t
        y.missing.t <- mu.t + sig.t * t(sd.mtx) %*% rnorm(ns, 0, 1)
        y.impute[missing.these, t] <- y.missing.t[missing.these]
        
      }
      
      # Only the sites/days with missing/thresholded observations are different from
      # the true y in y.imputed
      y <- y.impute

    }
    
    # update beta
    mmm <- rep(beta.m, p)
    vvv <- diag(p) / (beta.s^2)
    
    for (t in 1:nt) {
      taug.t <- sqrt(taug[, t])
      res.t  <- (y[, t] - lambda * zg[, t]) * taug.t
      x.t    <- x[, t, ] * taug.t
      ttt    <- t(x.t) %*% prec.cor
      vvv    <- vvv + ttt %*% x.t
      mmm    <- mmm + ttt %*% res.t
    }
    
    vvv  <- chol2inv(chol(vvv))
    mmm  <- vvv %*% mmm
    beta <- mmm + t(chol(vvv)) %*% rnorm(p)
    # print(mmm)
    # print(vvv)
    for (t in 1:nt) {
      x.beta[, t] <- x[, t, ] %*% beta
    }
    mu  <- x.beta + lambda * zg
    res <- y - mu
    
    # update partitions
    if (nknots > 1) {
      avgparts <- rep(0, nt)
      for (t in 1:nt) {
        att.w[1] <- att.w[1] + 1
        can.knots.con <- knots.con[, , t] + mh.w[1] * rnorm(2 * nknots)
        can.knots     <- pnorm(can.knots.con)
        can.knots[, 1] <- can.knots[, 1] * range + min(s[, 1])
        can.knots[, 2] <- can.knots[, 2] * range + min(s[, 2])
        can.g          <- mem(s, can.knots)
        can.taug       <- tau[can.g, t]
        can.zg         <- z[can.g, t]
        can.res        <- y - x.beta - lambda * can.zg
        
        if (temporalw & (t > 1)) {  # first day has mean 0: added for ts
          mean <- phi.w * knots.con[, , (t - 1)]
          sd   <- sqrt(1 - phi.w^2)
        } else {
          mean <- 0
          sd   <- 1
        }
        
        R <- -0.5 * quad.form(prec.cor, sqrt(can.taug) * can.res[, t]) +
              0.5 * quad.form(prec.cor, sqrt(taug[, t]) * res[, t]) +
              0.5 * sum(log(can.taug)) - 0.5 * sum(log(taug[, t])) +
              sum(dnorm(can.knots.con, mean, sd, log=T)) -  # added for ts
              sum(dnorm(knots.con[, , t], mean, sd, log=T))  # added for ts
        
        if (temporalw & (t < nt)) {  # the knot location on the next day is a part of the time series
          sd.next <- sqrt(1 - phi.w^2)
          knots.next <- knots.con[, , (t + 1)]
          R <- R + sum(dnorm(knots.next, (phi.w * can.knots.con), sd.next, log=T)) - 
                   sum(dnorm(knots.next, (phi.w * knots.con[, , t]), sd.next, log=T))
        }
        
        if (!is.na(R)) { if (log(runif(1)) < R) {
          knots.con[, , t] <- can.knots.con
          knots[, , t]     <- can.knots
          g[, t]           <- can.g
          taug[, t]        <- can.taug
          zg[, t]          <- can.zg
          acc.w[1]         <- acc.w[1] + 1
        }}
      }
    }  # fi nknots > 1
    
    # covariance parameters
    mu <- x.beta + lambda * zg
    res <- y - mu
    # update tau
    if (method == "gaussian") {  # single random effect for all days
      rss <- 0
      for (t in 1:nt) {
        res.t <- res[, t]
        rss <- rss + quad.form(prec.cor, sqrt(taug[, t]) * res.t)
      }
      tau <- rgamma(1, ns * nt / 2 + tau.alpha, rss / 2 + tau.beta)
      tau <- matrix(tau, nknots, nt)
      taug <- matrix(tau, ns, nt)
      
      # update tau.alpha and tau.beta
      a.star <- tau.beta.a + tau.alpha
      b.star <- tau.beta.b + tau[1, 1]
      tau.beta <- rgamma(1, a.star, b.star)
      
      lll <- mmm <- seq(0.5, 10, 0.1)
      for (l in 1:length(lll)) {
        lll[l] <- sum(dgamma(tau[1, 1], mmm[l], tau.beta, log=T))
      }
      tau.alpha <- sample(mmm, 1, prob=exp(lll - max(lll)))
    } else if (method == "t") { 
      if (nknots == 1) {
        for (t in 1:nt) {
          res.t <- res[, t]
          rss.t <- quad.form(prec.cor, res.t)
          
          aaa <- tau.alpha + 0.5 * ns
          bbb <- tau.beta + 0.5 * rss.t
          if (skew) {  # tau is also in z likelihood
            aaa <- aaa + 0.5
            bbb <- bbb + 0.5 * z[1, t]^2
          }
          
          if (!temporaltau) {  # conjugate
            tau[1, t] <- rgamma(1, aaa, bbb)
            taug[, t] <- tau[1, t]
          } else {  # not conjugate
            # TODO: time series update
          }  # fi temporaltau
        }  
      } else {  # nknots > 1
      	for (t in 1:nt) {
      	  res.t <- res[, t]
      	  cur.lly <- 0.5 * sum(log(taug[, t])) - 
      	             0.5 * quad.form(prec.cor, sqrt(taug[, t]) * res.t)
      	              
      	  for (k in 1:nknots) {
      	    these <- which(g[, t] == k)
      	    nparts <- length(these)
      	    nparts.tau[k, t] <- nparts
      	    
      	    if (nparts == 0) {
      	      aaa <- tau.alpha
      	      bbb <- tau.beta
      	      if (skew) {  # tau is also in z likelihood
      	        aaa <- aaa + 0.5
      	        bbb <- bbb + 0.5 * z[k, t]^2
      	      }
      	      
      	      if (!temporaltau) {
      	        tau[k, t] <- rgamma(1, aaa, bbb)
      	        if(tau[k, t] < 1e-6) {  # numerical stability
      	          tau[k, t] <- 1e-6
      	        }
      	      } else { # tau is not conjugate
      	        # TODO: time series update
      	      }
      	    } else {  # nparts > 0
      	      att.tau.ns[(nparts + 1)] <- att.tau.ns[(nparts + 1)] + 1
      	      att.tau[k, t] <- att.tau[k, t] + 1
      	      
      	      aaa <- tau.alpha + 0.5 * nparts
      	      bbb <- tau.beta + 0.5 * quad.form(prec.cor[these, these], res.t[these])
      	      
      	      if (skew) {  # tau is also in z likelihood
      	        aaa <- aaa + 0.5
      	        bbb <- bbb + 0.5 * z[k, t]^2
      	      }
      	      
      	      # this posterior is conjugate when nparts -> ns 
      	      # when nparts -> 1, want a wider candidate
      	      aaa <- aaa / mh.tau.ns[nparts + 1]
      	      bbb <- bbb / mh.tau.ns[nparts + 1]
      	      if(is.nan(bbb)) {
      	        print(tau)
      	        print(z)
      	        print(aaa)
      	        print(y[1:5, ])
      	        print(res[1:5, ])
      	        print(mu[1:5, ])
      	        print(beta)
      	        print(quad.form(prec.cor[these, these], res.t[these]))
      	      }
      	      can.tau <- tau[, t]
      	      can.tau[k] <- rgamma(1, aaa, bbb)
      	      if (can.tau[k] < 1e-6) {
      	        can.tau[k] <- 1e-6
      	      }
      	      
      	      can.taug <- can.tau[g[, t]]
      	      
      	      can.lly <- 0.5 * sum(log(can.taug)) - 
      	                 0.5 * quad.form(prec.cor, sqrt(can.taug) * res.t)
      	      
      	      if (skew) {
      	        cur.llz <- 0.5 * log(tau[k, t]) - 0.5 * tau[k, t] * z[k, t]^2
      	        can.llz <- 0.5 * log(can.tau[k]) - 0.5 * can.tau[k] * z[k, t]^2
      	      } else {
      	        cur.llz <- can.llz <- 0
      	      }
      	      
      	      R <- can.lly - cur.lly + can.llz - cur.llz +
      	           tryCatch({  # candidate is non-symmetric
                     dgamma(tau[k, t], aaa, bbb, log=TRUE)}, 
                     warning = function(w) {
                     print(paste("knot", k, ", day", t))
                     print(paste("aaa =", aaa))
                     print(paste("bbb =", bbb))
                     print(paste("tau[k, t] =", tau[k, t]))
                     print(paste("g[, t] =", g[, t]))
                     print(paste("mh.tau.ns[mh.idx] =", mh.tau.ns[mh.idx]))
                     print(paste("tau.beta =", tau.beta))
                   }) -
                   dgamma(can.tau[k], aaa, bbb, log=TRUE)

      	      if (!temporaltau) {
      	        R <- R + dgamma(can.tau[k], tau.alpha, tau.beta, log=T) -
      	                 dgamma(tau[k, t], tau.alpha, tau.beta, log=T)
      	      } else {  # prior changes
      	        # TODO: time series update
      	      }
      	      
      	      if (!is.na(R)) { if (log(runif(1)) < R) {
      	        acc.tau.ns[(nparts + 1)] <- acc.tau.ns[(nparts + 1)] + 1
      	        acc.tau[k, t] <- acc.tau[k, t] + 1
      	        tau[, t] <- can.tau
      	        taug[, t] <- can.taug
      	        cur.lly <- can.lly
      	      }}
      	      
      	    }  # fi nparts
      	  }  # end k
      	}  # end t
      }  # fi nknots > 1
      
      # update hyperparameters
      if (fixhyper) {
        tau.alpha <- 3
        tau.beta  <- 8
      } else {
        a.star <- tau.beta.a + tau.alpha * nknots * nt
        b.star <- tau.beta.b + sum(tau)
        tau.beta <- rgamma(1, a.star, b.star)
        
        lll <- mmm <- seq(0.1, 10, 0.1)
        for (l in 1:length(lll)) {
          lll[l] <- sum(dgamma(tau, mmm[l], tau.beta, log=T))
        }
        tau.alpha <- sample(mmm, 1, prob=exp(lll - max(lll)))
      }
    }  # fi method == t
    
    mu  <- x.beta + lambda * zg
    res <- y - mu
    # update rho and nu and gamma
    if (rho.prior == "disc") {  # only update rho and alpha
      lll.rho <- rep(NA, length(rhos))
      
      # storage for possible covariance parts
      can.prec.cor <- array(NA, dim=c(ns, ns, length(rhos)))
      can.logdet.prec <- rep(NA, length(rhos))
      can.rss <- matrix(NA, nrow=length(rhos), ncol=nt)
      for (l in 1:length(rhos)) {
      	D <- 1 - gamma + gamma * C.values[, l]  # eigenvalues with gamma
      	can.prec.cor[, , l] <- quad.tform(diag(1 / D), C.vectors[, , l])
      	can.logdet.prec[l] <- -0.5 * sum(log(D))
      	for (t in 1:nt) {
      	  can.rss[l, t] <- quad.form(can.prec.cor[, , l], sqrt(taug[, t]) * res[, t])
      	}
      	lll.rho[l] <- nt * can.logdet.prec[l] - 0.5 * sum(can.rss[l, ])
      }
      rho.idx <- sample(length(rhos), 1, prob=exp(lll.rho - max(lll.rho)))
      rho <- rhos[rho.idx]
      
      # update cov matrix
      C <- exp(-(d / rho))
      prec.cor <- can.prec.cor[, , rho.idx]
      logdet.prec <- can.logdet.prec[rho.idx]
      cur.rss <- can.rss[rho.idx, ]
            
      # update gamma
      att.gamma <- att.gamma + 1
    
      norm.gamma <- qnorm(gamma)
      can.norm.gamma <- rnorm(1, norm.gamma, mh.gamma)
      can.gamma <- pnorm(can.norm.gamma)
    
      can.C <- CorFx(d=d, gamma=can.gamma, rho=rho, nu=nu)
      can.CC <- tryCatch(chol.inv(can.C, inv=T, logdet=T),
                         error = function(e) {
                           eig.inv(can.C, inv=T, logdet=T, mtx.sqrt=T)
                         })
      can.sd.mtx <- can.CC$sd.mtx
      can.prec.cor <- can.CC$prec
      can.logdet.prec <- can.CC$logdet.prec  # this is the sqrt of logdet.prec
    
      can.rss <- rep(NA, nt)
      for (t in 1:nt) {
        can.rss[t] <- quad.form(can.prec.cor, sqrt(taug[, t]) * res[, t])
        cur.rss[t] <- quad.form(can.prec.cor, sqrt(taug[, t]) * res[, t])
      }
    
      R <- -0.5 * sum(can.rss - cur.rss) + 
            nt * (can.logdet.prec - logdet.prec) +
            dnorm(can.norm.gamma, mean=gamma.m, sd=gamma.s, log=T) - 
            dnorm(norm.gamma, mean=gamma.m, sd=gamma.s, log=T)
    
      if (!is.na(R)) { if (log(runif(1)) < R) {
        gamma <- can.gamma
        C <- can.C
        sd.mtx <- can.sd.mtx
        prec.cor <- can.prec.cor
        logdet.prec <- can.logdet.prec
        cur.rss <- can.rss
        acc.gamma <- acc.gamma + 1
      }}
    
      if ((att.gamma > 50) & (iter < (burn / 2))) {
        if (acc.gamma / att.gamma < 0.25) { mh.gamma <- mh.gamma * 0.8 }
        if (acc.gamma / att.gamma > 0.50) { mh.gamma <- mh.gamma * 1.2 }
        acc.gamma <- att.gamma <- 0
      }

    } else {  # using a truncated normal candidate
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
    
      can.C <- CorFx(d=d, gamma=gamma, rho=can.rho, nu=can.nu)
      can.CC <- tryCatch(chol.inv(can.C, inv=T, logdet=T),
                         error = function(e) {
                           tryCatch(eig.inv(can.C, inv=T, logdet=T, mtx.sqrt=T),
                                    error = function(e) {
                                      print(paste("can.rho =", can.rho))
                                      print(paste("can.nu =", can.nu))
                                    })
                         })
      can.sd.mtx <- can.CC$sd.mtx
      can.prec.cor <- can.CC$prec
      can.logdet.prec <- can.CC$logdet.prec  # this is the sqrt of logdet.prec
    
      can.rss <- rep(NA, nt)
      cur.rss <- rep(NA, nt)
      for (t in 1:nt) {
        can.rss[t] <- quad.form(can.prec.cor, sqrt(taug[, t]) * res[, t]) 
        cur.rss[t] <- quad.form(prec.cor, sqrt(taug[, t]) * res[, t])
      }
    
      R <- -0.5 * sum(can.rss - cur.rss) + 
            nt * (can.logdet.prec - logdet.prec) + 
            dnorm(can.logrho, logrho.m, logrho.s, log=T) - 
            dnorm(logrho, logrho.m, logrho.s, log=T) + 
            dnorm(can.lognu, lognu.m, lognu.s, log=T) - 
            dnorm(lognu, lognu.m, lognu.s, log=T)
      
      if (upper.logrho < 1) {  # candidate is not symmetric
        R <- R + dnorm(logrho, logrho, mh.rho, log=T) - 
                 dnorm(can.logrho, logrho, mh.rho, log=T)
      }
      if (upper.lognu < 1) {  # candidate is not symmetric
        R <- R + dnorm(lognu, lognu, mh.nu, log=T) - 
                 dnorm(can.lognu, lognu, mh.nu, log=T)
      }
    
      if (!is.na(R)) { if (log(runif(1)) < R) {
        rho <- can.rho
        nu <- can.nu
        C <- can.C
        sd.mtx <- can.sd.mtx
        prec.cor <- can.prec.cor
        logdet.prec <- can.logdet.prec
        cur.rss <- can.rss
        acc.rho <- acc.rho + 1
        acc.nu  <- acc.nu + 1
      }}
           
      # gamma
      att.gamma <- att.gamma + 1
      
      norm.gamma <- qnorm(gamma)
      can.norm.gamma <- rnorm(1, norm.gamma, mh.gamma)
      can.gamma <- pnorm(can.norm.gamma)
    
      can.C <- CorFx(d=d, gamma=can.gamma, rho=rho, nu=nu)
      can.CC <- tryCatch(chol.inv(can.C, inv=T, logdet=T),
                         error = function(e) {
                           tryCatch(eig.inv(can.C, inv=T, logdet=T, mtx.sqrt=T),
                                    error = function(e) {
                                      print(paste("can.gamma =", can.gamma))
                                    })
                         })
      can.sd.mtx <- can.CC$sd.mtx
      can.prec.cor <- can.CC$prec
      can.logdet.prec <- can.CC$logdet.prec  # this is the sqrt of logdet.prec
    
      can.rss <- rep(NA, nt)
      cur.rss <- rep(NA, nt)
      for (t in 1:nt) {
        can.rss[t] <- quad.form(can.prec.cor, sqrt(taug[, t]) * res[, t]) 
        cur.rss[t] <- quad.form(prec.cor, sqrt(taug[, t]) * res[, t])
      }
    
      R <- -0.5 * sum(can.rss - cur.rss) + 
            nt * (can.logdet.prec - logdet.prec) + 
            dnorm(can.norm.gamma, mean=gamma.m, sd=gamma.s, log=T) - 
            dnorm(norm.gamma, mean=gamma.m, sd=gamma.s, log=T)
    
      if (!is.na(R)) { if (log(runif(1)) < R) {
        acc.gamma <- acc.gamma + 1
        gamma <- can.gamma
        C <- can.C
        sd.mtx <- can.sd.mtx
        prec.cor <- can.prec.cor
        logdet.prec <- can.logdet.prec
        cur.rss <- can.rss
      }}
            
      if ((att.rho > 50) & (iter < (burn / 2))) {
        if (acc.rho / att.rho < 0.25) { mh.rho <- mh.rho * 0.8 }
        if (acc.rho / att.rho > 0.50) { mh.rho <- mh.rho * 1.2 }
        acc.rho <- att.rho <- 0
      }
    
      if ((att.nu > 50) & (iter < (burn / 2))) {
        if (acc.nu / att.nu < 0.25) { mh.nu <- mh.nu * 0.8 }
        if (acc.nu / att.nu > 0.50) { mh.nu <- mh.nu * 1.2 }
        acc.nu <- att.nu <- 0
      }
    
      if ((att.gamma > 50) & (iter < (burn / 2))) {
        if (acc.gamma / att.gamma < 0.25) { mh.gamma <- mh.gamma * 0.8 }
        if (acc.gamma / att.gamma > 0.50) { mh.gamma <- mh.gamma * 1.2 }
        acc.gamma <- att.gamma <- 0
      }
      
      if (fixhyper) {
        nu <- 0.5
        rho <- 1
        gamma <- 0.9
      }

    }
    
    # update skew parameters: lambda and z
    if (skew) {
      # lambda
      mmm <- lambda.m
      vvv <- 1 / (lambda.s^2)
      
      for (t in 1:nt) {
        taug.t <- sqrt(taug[, t])
        res.t  <- (y[, t] - x.beta[, t]) * taug.t
        z.t    <- zg[, t] * taug.t
        ttt    <- z.t %*% prec.cor
        vvv    <- vvv + ttt %*% z.t
        mmm    <- mmm + ttt %*% res.t               
      }
      
      vvv <- 1 / vvv
      mmm <- vvv * mmm
      lambda <- rnorm(1, mmm, sqrt(vvv))
      
      mu <- x.beta + lambda * zg
      res <- y - mu
      
      # z
      if (!temporalz) {
        if (nknots == 1) {
          for (t in 1:nt) {
          	res.t <- (y[, t] - x.beta[, t])
            mmm <- lambda * tau[1, t] * sum(prec.cor %*% res.t)
            vvv <- tau[1, t] + lambda^2 * tau[1, t] * sum(prec.cor)
            
            vvv <- 1 / vvv
            mmm <- vvv * mmm
            z[1, t] <- abs(rnorm(1, mmm, sqrt(vvv)))
            zg[, t] <- z[1, t]
          }
        } else {  # nknots > 1
          for (t in 1:nt) {
          	taug.t <- sqrt(taug[, t])
          	for (k in 1:nknots) {
              these <- which(g[, t] == k)
              r.1 <- (y[these, t] - x.beta[these, t]) * taug.t[these]
              r.2 <- (y[-these, t] - mu[-these, t]) * taug.t[-these]
              
              prec.11 <- prec.cor[these, these, drop=F]
              prec.21 <- prec.cor[-these, these, drop=F]
              
              mmm <- lambda * sqrt(tau[k, t]) * sum(r.1 %*% prec.11 + r.2 %*% prec.21)
              vvv <- tau[k, t] + lambda^2 * tau[k, t] * sum(prec.11)
              
              vvv <- 1 / vvv
              mmm <- vvv * mmm
              z[k, t] <- abs(rnorm(1, mmm, sqrt(vvv)))
              zg[these, t] <- z[k, t]
              mu[, t] <- x.beta[, t] + lambda * zg[, t]
            }
          }
        }  # fi nknots > 1
      } else {
        # TODO: time series z
      }  # fi temporalz
      mu <- x.beta + lambda * zg
      res <- y - mu
    }  # fi skew
        
  }  # end nthin
  
  mu <- x.beta + lambda * zg
  res <- y - mu
  # predictions
  if (predictions) {
  	if (cov.model == "matern") {
  	  s.11 <- gamma * simple.cov.sp(D=d11, sp.type="matern", sp.par=c(1, rho), error.var=0, 
                                    smoothness=nu, finescale.var=0)
      s.12 <- gamma * simple.cov.sp(D=d12, sp.type="matern", sp.par=c(1, rho), error.var=0,
                                    smoothness=nu, finescale.var=0)
    } else {
      s.11 <- gamma * matrix(exp(-d11 / rho), np, np)
      s.12 <- gamma * matrix(exp(-d12 / rho), np, ns)
    }
    diag(s.11) <- 1
    s.12.22.inv <- s.12 %*% prec.cor
    corp <- s.11 - s.12.22.inv %*% t(s.12)
    corp.sd.mtx <- tryCatch(chol(corp),  # only want the cholesky factor
                            error = function(e) {
                              eig.inv(corp, inv=F, logdet=F, mtx.sqrt=T)$sd.mtx
                            })
    
    yp <- matrix(NA, np, nt)

    for (t in 1:nt) {
      xp.beta  <- x.pred[, t, ] %*% beta
      if (nknots == 1) {
        gp <- rep(1, np)
      } else {
        gp <- mem(s.pred, knots[, , t])  # find the right partition
      }
      zgp    <- z[gp, t]
      siggp  <- 1 / sqrt(tau[gp, t])  # get the partition's standard deviation
      taug.t <- sqrt(taug[, t])
      mup <- xp.beta + lambda * zgp + siggp * s.12.22.inv %*% (taug.t * res[, t])
      
      yp[, t] <- mup + siggp * t(corp.sd.mtx) %*% rnorm(np, 0, 1)
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
  
  if (predictions) {
    y.pred[iter, , ] <- yp
  }
  if (skew) {
    keepers.lambda[iter] <- lambda
    keepers.z[iter, , ]   <- z
  }
  if (nknots > 1) {
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
 	}
  
  	acc.rate.rho <- round(acc.rho / att.rho, 3)
  	acc.rate.nu <- round(acc.nu / att.nu, 3)
  	acc.rate.tau <- round(acc.tau / att.tau, 3)
  	acc.rate.gamma <- round(acc.gamma / att.gamma, 3)
  	
  	if (iter < burn) {
  	  begin <- max(1, (iter - 2000))
  	} else {
  	  begin <- burn
  	}
  	
    if (iterplot) {
      if (skew) {
        par(mfrow=c(3, 6))
      } else {
        par(mfrow=c(3, 6))
      }
      
      plot(keepers.beta[begin:iter, 1], type="l")
      if (p > 1) {
        plot(keepers.beta[begin:iter, 2], type="l")
      } else if (p > 2) {
        plot(keepers.beta[begin:iter, 3], type="l")
      }
      if (temporalw) {
        title.phi.w <- paste("acc =", acc.rate.phi.w)
        plot(keepers.phi.w[begin:iter], type="l", main=title.phi.w)
      }
      if (temporalz) {
        title.phi.z <- paste("acc =", acc.rate.phi.z)
      	plot(keepers.phi.z[begin:iter], type="l", main=title.phi.z)
      }
      if (temporaltau) {
        title.phi.tau <- paste("acc =", acc.rate.phi.tau)
      	plot(keepers.phi.tau[begin:iter], type="l", main=title.phi.tau)
      }
      
      plot(keepers.tau.alpha[begin:iter], type="l")
      plot(keepers.tau.beta[begin:iter], type="l")
      
      title.rho <- paste("acc =", acc.rate.rho)
      if (mh.rho < 0.00001) {
      	xlab.rho <- "<0.00001"
      } else if (mh.rho < 10000) {
      	xlab.rho <- round(mh.rho, 5)
      } else {
      	xlab.rho <- "> 10000"
      }
      plot(keepers.rho[begin:iter], type="l", main=title.rho, xlab=xlab.rho)

      title.nu <- paste("acc =", acc.rate.nu)
      if (mh.nu < 0.00001) {
      	xlab.nu <- "<0.00001"
      } else if (mh.nu < 10000) {
      	xlab.nu <- round(mh.nu, 5)
      } else {
      	xlab.nu <- "> 10000"
      }
      plot(keepers.nu[begin:iter], type="l", main=title.nu, xlab=xlab.nu)
      
      title.gamma <- paste("acc =", acc.rate.gamma)
      plot(keepers.gamma[begin:iter], type="l", main=title.gamma)
      
      if (skew) {
        plot(keepers.lambda[begin:iter], type="l")
        
        if (temporalz) {
          title.z.1 <- paste("acc =", acc.rate.z[1, 1])
          title.z.2 <- paste("acc =", acc.rate.z[1, 10])
          title.z.3 <- paste("acc =", acc.rate.z[1, 21])
        } else {
          title.z.1 <- ""
          title.z.2 <- ""
          title.z.3 <- ""
        }
        plot(keepers.z[begin:iter, 1, 1], type="l", main=title.z.1)
        plot(keepers.z[begin:iter, 1, 10], type="l", main=title.z.2)
        plot(keepers.z[begin:iter, 1, 21], type="l", main=title.z.3)
      } 
      
      if (temporaltau) {
        nparts.1 <- nparts.2 <- nparts.3 <- 1
      } else {
        nparts.1 <- length(which(g[, 1] == 1))
        nparts.2 <- length(which(g[, 10] == 1))
        nparts.3 <- length(which(g[, 21] == 1))
      }      
      mh.disp.1 <- round(mh.tau.ns[(nparts.1 + 1)], 2)
      mh.disp.2 <- round(mh.tau.ns[(nparts.2 + 1)], 2)
      mh.disp.3 <- round(mh.tau.ns[(nparts.3 + 1)], 2)
      title.tau.1 <- paste("acc = ", acc.rate.tau[1, 1])
      title.tau.2 <- paste("acc = ", acc.rate.tau[1, 10])
      title.tau.3 <- paste("acc = ", acc.rate.tau[1, 21])
      plot(keepers.tau[begin:iter, 1, 1], type="l", main=title.tau.1, 
           ylab="tau 1,1", xlab=paste(nparts.1, ", ", mh.disp.1))
      plot(keepers.tau[begin:iter, 1, 10], type="l", main=title.tau.2, 
           ylab="tau 1, 10", xlab=paste(nparts.2, ", ", mh.disp.2))
      plot(keepers.tau[begin:iter, 1, 21], type="l", main=title.tau.3, 
           ylab="tau 1, 21", xlab=paste(nparts.3, ", ", mh.disp.3))
     
      if (nknots > 1) {
      	if (temporaltau) {
      	  nparts.4 <- nparts.5 <- nparts.6 <- 1
      	} else {
      	  nparts.4 <- length(which(g[, 1] == 2))
      	  nparts.5 <- length(which(g[, 10] == 2))
      	  nparts.6 <- length(which(g[, 21] == 2))
      	}
        
        mh.disp.4 <- round(mh.tau.ns[(nparts.4 + 1)], 2)
        mh.disp.5 <- round(mh.tau.ns[(nparts.5 + 1)], 2)
        mh.disp.6 <- round(mh.tau.ns[(nparts.6 + 1)], 2)
        title.tau.4 <- paste("acc = ", acc.rate.tau[2, 1])
        title.tau.5 <- paste("acc = ", acc.rate.tau[2, 10])
        title.tau.6 <- paste("acc = ", acc.rate.tau[2, 21])
        plot(keepers.tau[begin:iter, 2, 1], type="l", main=title.tau.4, 
             ylab="tau 2, 1", xlab=paste(nparts.4, ", ", mh.disp.4))
        plot(keepers.tau[begin:iter, 2, 10], type="l", main=title.tau.5, 
             ylab="tau 2, 10", xlab=paste(nparts.5, ", ", mh.disp.5))
        plot(keepers.tau[begin:iter, 2, 21], type="l", main=title.tau.6, 
             ylab="tau 2, 21", xlab=paste(nparts.6, ", ", mh.disp.6))
      }    
      
    }
    
    cat("\t iter", iter, "\n")
  }
  
  
  } #end iters

if (nknots == 1) {
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
} else {
  y.pred <- y.pred[return.iters, , ]
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

