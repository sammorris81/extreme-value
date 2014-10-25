#########################################################################
# MCMC 
#
# TODO: Add in model description here
#
#ASSUMES x and y are in [0,1]^2
#########################################################################

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL, 
                 thresh.all=0, thresh.quant=T, nknots=1, keep.knots=F,
                 iters=5000, burn=1000, update=100, thin=1,
                 iterplot=F, plotname=NULL, method="t", 
                 # just to debug temporal parts.
                 temporalw=F, temporaltau=F, temporalz=F,  # eventually change to temporal=F
                 # initial values
                 beta.init=NULL, 
                 tau.init=2,
                 tau.alpha.init=0.1, 
                 tau.beta.init=0.1,
                 rho.init=0.5, 
                 nu.init=0.5,
                 alpha.init=0.5,
                 # priors
                 beta.m=0, beta.s=10, 
                 tau.alpha.m=0, tau.alpha.s=1, 
                 tau.beta.a=0.1, tau.beta.b=0.1,
                 logrho.m=0, logrho.s=10,
                 lognu.m=-1.2, lognu.s=1,
                 alpha.m=0, alpha.s=1,
                 # covariance model
                 cov.model="matern",  # or "exponential"
                 rho.prior="cont",  # or "disc"
                 # skew inits
                 z.init=1, z.alpha.init=0,
                 # skew priors
                 z.alpha.m=0, z.alpha.s=2, skew=T,
                 thresh.site.specific=F, thresh.site=NULL
        ){
    
  library(geoR)
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
  
  zg <- matrix(0, ns, nt)
  if (length(z.init) == 1 && skew) {
    cat("\t initializing all z terms to", z.init, "\n")
  }
  z <- matrix(z.init, nknots, nt)
  for (t in 1:nt) {
    zg[, t] <- z[g[, t], t]
  }
  
  if (skew) {
    z.alpha <- z.alpha.init
  } else {
  	if (z.alpha.init != 0) {
  	  warning("z.alpha.init being ignored since skew=F")
  	}
    z.alpha <- 0
  }
    
  # easier to keep calculations in the precision scale for MCMC
  sigma2    <- 1 / tau
  sigma2g   <- 1 / taug
  tau.alpha <- tau.alpha.init
  tau.beta  <- tau.beta.init
  
  # initialize spatial covariance
  rho    <- rho.init
  logrho <- log(rho)
  nu     <- nu.init
  lognu  <- log(nu)
  alpha  <- alpha.init
  fixnu <- F
  
  if (cov.model == "exponential") {
    nu <- 0.5
    fixnu <- T
  }
  
  if (rho.prior == "cont") {
    C <- CorFx(d=d, alpha=alpha, rho=rho, nu=nu)
    CC <- tryCatch(chol.inv(C, inv=T, logdet=T),
                   error = function(e) {
                     eig.inv(C, inv=T, logdet=T, mtx.sqrt=T)
                   })
    chol.C <- CC$sd.mtx
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
   	D <- 1 - alpha + alpha * C.values[, rho.idx]  # eigenvalues with alpha
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
    mh.tau.ns <- rep(0.1, 10)
  } else {
  	mh.tau.ns <- rep(1, 10)
  }
  if (temporalz) {
  	z.star <- z  # need a place to keep track of normal values for time series
    phi.z <- 0.5
    acc.z <- att.z <- mh.z <- matrix(1, nknots, nt)
    acc.phi.z <- att.phi.z <- mh.phi.z <- 1
  }
    
  # MH tuning params
  acc.w      <- att.w      <- mh.w     <- rep(0.1, nt)  # knot locations
  acc.delta  <- att.delta  <- mh.delta <- 0.1  
  acc.rho    <- att.rho    <- mh.rho   <- 0.1
  acc.nu     <- att.nu     <- mh.nu    <- 0.1
  acc.alpha  <- att.alpha  <- mh.alpha <- 0.5
  
  # candidate distributions for tau depend on the percentage of sites in the 
  # partition. 
  acc.tau.ns <- att.tau.ns <- rep(1, ns)
  avg.sites.part <- 1 / nknots
  mh.tau.parts <- c(0, (avg.sites.part / 4), (avg.sites.part / 2), avg.sites.part, 
                    (2 * avg.sites.part), (3 * avg.sites.part), (4 * avg.sites.part))
  acc.tau    <- att.tau   <- matrix(1, nrow=nknots, ncol=nt) 
  mh.tau <- matrix(0.05, nknots, nt)
  nparts.tau <- matrix(1, nrow=nknots, ncol=nt)
       
  # storage
  keepers.tau       <- array(NA, dim=c(iters, nknots, nt))
  keepers.beta      <- matrix(NA, nrow=iters, ncol=p)
  keepers.ll        <- matrix(NA, nrow=iters, ncol=nt)
  keepers.tau.alpha <- rep(NA, iters)
  keepers.tau.beta  <- rep(NA, iters)
  keepers.nparts    <- matrix(NA, nrow=iters, ncol=nt)
  keepers.rho       <- rep(NA, iters)
  keepers.nu        <- rep(NA, iters)
  keepers.alpha     <- rep(NA, iters)
  keepers.z         <- array(NA, dim=c(iters, nknots, nt))
  keepers.z.alpha   <- rep(NA, iters)
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
    
    # impute data below threshold
    if (thresh.all != 0) {
      mu <- x.beta + z.alpha * zg
      thresh.mtx.fudge <- 0.99999 * thresh.mtx  # numerical stability
      y.imputed <- matrix(y, ns, nt)
      
      if (nu == 0.5) {  # quicker than matern function
        cor <- alpha * exp(-(d / rho))
      } else {
        cor <- alpha * matern(u=d, phi=rho, kappa=nu)  # only for the spatial error
      }
      for (t in 1:nt) {
        these.thresh.obs <- thresh.obs[, t]
        these.missing.obs <- missing.obs[, t]

        # spatial error
        sig.t <- 1 / sqrt(taug[, t])
        sd.mtx <- tryCatch(chol(cor),  # only want the cholesky factor
                           error = function(e) {
                             eig.inv(cor, inv=F, logdet=F, mtx.sqrt=T)$sd.mtx
                           })
        theta.t <- sig.t * t(sd.mtx) %*% rnorm(ns, 0, 1)  # generate for all sites
        
        # nugget error
        # new expected value and standard deviation
        e.y <- mu[, t] + theta.t
        s.y <- sqrt(1 - alpha) * sig.t
        upper.y <- thresh.mtx[, t]
        
        y.impute.t <- tryCatch(rTNorm(mn=e.y, sd=s.y, lower=-Inf, upper=upper.y),
                               warning = function(e) {
                               	cat("sig.t = ", sig.t, "\n")
                               	cat("e.y = ", e.y, "\n")
                               	cat("s.y = ", s.y, "\n")
                               	cat("alpha = ", alpha, "\n")
                               })
       
        # if any y.impute.t come back as -Inf, it's because P(Y < T) = 0
        usefudge <- y.impute.t == -Inf
        y.impute.t[usefudge] <- thresh.mtx.fudge[usefudge, t]
        y.imputed[these.thresh.obs, t] <- y.impute.t[these.thresh.obs]
       
        # we already know the spatial part through the theta term
        y.missing.t <- rnorm(n=ns, mean=e.y, sd=s.y)
        y.imputed[these.missing.obs, t] <- y.missing.t[these.missing.obs]
      }
      
      # Only the sites/days with missing/thresholded observations are different from
      # the true y in y.imputed
      y <- y.imputed

    }
    
    # update beta
    vvv <- diag(p) / beta.s^2
    mmm <- rep(beta.m, p)  
    for (t in 1:nt) {
       x.t    <- x[, t, ]
       taug.t <- sqrt(taug[, t])
       prec.t <- sweep(prec.cor, 2, taug.t, "*") * taug.t
       ttt    <- t(x.t) %*% prec.t
       vvv    <- vvv + ttt %*% x.t
       mmm    <- mmm + ttt %*% (y[, t] - z.alpha * zg[, t])
    }
    
    vvv <- chol2inv(chol(vvv))
    beta <- vvv %*% mmm + t(chol(vvv)) %*% rnorm(p)
      
    for (t in 1:nt) {
      x.beta[, t] <- x[, t, ] %*% beta
    }
    mu  <- x.beta + z.alpha * zg
    
    # update partitions
    res <- y - mu
    if (nknots > 1) {
      avgparts <- rep(0, nt)
      for (t in 1:nt) {
        att.w[1]       <- att.w[1] + 1
        can.knots.con  <- knots.con[, , t] + mh.w[1] * rnorm(2 * nknots)
        can.knots      <- pnorm(can.knots.con)
        can.knots[, 1] <- can.knots[, 1] * range + min(s[, 1])
        can.knots[, 2] <- can.knots[, 2] * range + min(s[, 2])
        can.g          <- mem(s, can.knots)
        can.taug       <- tau[can.g, t]
        can.zg         <- z[can.g, t]
        can.res        <- y - x.beta - z.alpha * can.zg   
        
        if (temporalw & (t > 1)) {  # first day has mean 0: added for ts
          mean <- phi.w * knots.con[, , (t-1)]
          sd   <- sqrt(1 - phi.w^2)
        } else {
          mean <- 0
          sd   <- 1
        }

        R <- -0.5 * quad.form(prec.cor, sqrt(can.taug) * can.res[, t]) +
              0.5 * quad.form(prec.cor, sqrt(taug[, t]) * res[, t]) +
              0.5 * sum(log(can.taug)) -
              0.5 * sum(log(taug[, t])) +
              sum(dnorm(can.knots.con, mean, sd, log=T)) -  # edited for ts
              sum(dnorm(knots.con[, , t], mean, sd, log=T))  # edited for ts
        
        if (temporalw & (t < nt)) {
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
        
        # find the avg num of partitions per day
        for (k in 1:nknots) {
          if (sum(g[, t] == k) != 0) {
            avgparts[t] <- avgparts[t] + 1
          } 
        }
      }  # end t in 1:nt
      
      if (temporalw) {
        att.phi.w <- att.phi.w + 1
        phi.con.w <- qnorm((phi.w + 1) / 2)  # transform to R
        can.phi.con.w <- rnorm(1, phi.con.w, mh.phi.w)  # draw candidate
        can.phi.w <- 2 * pnorm(can.phi.con.w) - 1  # transform back to (-1, 1)
        
        can.ll <- cur.ll <- 0
        for (t in 2:nt) {  # no phi.w for day 1
          cur.mean <- phi.w * knots.con[, , (t - 1)]
          can.mean <- can.phi.w * knots.con[, , (t - 1)]
          cur.sd <- sqrt(1 - phi.w^2)
          can.sd <- sqrt(1 - can.phi.w^2)
          can.ll <- can.ll + sum(dnorm(knots.con[, , t], can.mean, can.sd, log=T))
          cur.ll <- cur.ll + sum(dnorm(knots.con[, , t], cur.mean, cur.sd, log=T))
        }
        
        R <- can.ll - cur.ll +
             dnorm(can.phi.con.w, log=T) - dnorm(phi.con.w, log=T)
             
        if (!is.na(R)) { if (log(runif(1)) < R) {
          acc.phi.w <- acc.phi.w + 1
          phi.w <- can.phi.w
        } }
        
        if ((att.phi.w > 50) & (iter < (burn / 2))) {
          if (acc.phi.w / att.phi.w < 0.25) { mh.phi.w <- mh.phi.w * 0.8 }
          if (acc.phi.w / att.phi.w > 0.50) { mh.phi.w <- mh.phi.w * 1.2 }
          acc.phi.w <- att.phi.w <- 0
        }
        
      }  # fi temporal
    }  # fi nknots > 1
    
    #### Spatial correlation
    mu <- x.beta + z.alpha * zg
    res <- y - mu

    # update tau
    curll <- rep(0,nt)
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
    } else if (method == "t"){
      if (temporaltau) {
        ts.tau <- ts.sample.tau(tau=tau, acc.tau=acc.tau, att.tau=att.tau, mh.tau=mh.tau, 
                                acc.tau.ns=acc.tau.ns, att.tau.ns=att.tau.ns, 
                                mh.tau.ns=mh.tau.ns, mh.tau.parts=mh.tau.parts,
                                taug=taug, phi=phi.tau, att.phi=att.phi.tau, 
                                acc.phi=acc.phi.tau, mh.phi=mh.phi.tau,  
                                tau.alpha=tau.alpha, tau.beta=tau.beta,
                                res=res, prec.cor=prec.cor, g=g, z=z)

        tau <- ts.tau$tau
        att.tau <- ts.tau$att.tau
        acc.tau <- ts.tau$acc.tau
        phi.tau <- ts.tau$phi
        att.phi.tau <- ts.tau$att.phi
        acc.phi.tau <- ts.tau$acc.phi
        if ((att.phi.tau > 50) & (iter < (burn / 2))) {
          if (acc.phi.tau / att.phi.tau < 0.25) { mh.phi.tau <- mh.phi.tau * 0.8 }
          if (acc.phi.tau / att.phi.tau > 0.50) { mh.phi.tau <- mh.phi.tau * 1.2 }
          acc.phi.tau <- att.phi.tau <- 0
        }

      } else {
        if (nknots == 1) {
          for (t in 1:nt) {
            res.t <- res[, t]
            rss.t <- quad.form(prec.cor, res.t)
            tau[1, t] <- rgamma(1, tau.alpha + 0.5 * ns, tau.beta + 0.5 * rss.t)
            taug[, t] <- tau[1, t]
          }  # end t in 1:nt

        } else {

          for(t in 1:nt){
            res.t <- res[, t]
            curll[t] <- 0.5 * sum(log(taug[, t])) -
                        0.5 * quad.form(prec.cor, sqrt(taug[, t]) * res.t)

            for (k in 1:nknots) {
              these  <- which(g[, t] == k)
              nparts <- length(these)
              # if (k == 1 & t == 10) {
                # print("knot 1 day 10")
                # print(paste("nparts is", nparts))
              # }
              nparts.tau[k, t] <- nparts
              
              if(nparts==0){
                tau[k, t] <- rgamma(1, tau.alpha, tau.beta)
                if (tau[k, t] < 1e-6) {
                  tau[k, t] <- 1e-6
                }
              } else {      
                mh.idx <- get.tau.mh.idx(nparts, ns, mh.tau.parts)
                att.tau.ns[mh.idx] <- att.tau.ns[mh.idx] + 1
                att.tau[k, t] <- att.tau[k, t] + 1
                aaa <- nparts / 2 + tau.alpha
                bbb <- quad.form(prec.cor[these, these], res.t[these]) / 2 + tau.beta
                # if ((k == 1 & t == 1) | (k == 1 & t == 10) | (k == 1 & t == 21)) {
                  # print(paste("MH.idx =", mh.idx))
                  # print(paste("MH s =", mh.tau.ns[mh.idx]))
                  # print(paste("nparts =", nparts))
                # }

                aaa <- aaa / mh.tau.ns[mh.idx]
                bbb <- bbb / mh.tau.ns[mh.idx]
                # aaa <- aaa / mh.tau[k, t]
                # bbb <- bbb / mh.tau[k, t]

                cantau    <- tau[, t]
                cantau[k] <- rgamma(1, aaa, bbb)
                if (cantau[k] < 1e-6) {
                  cantau[k] <- 1e-6
                }
                cantaug   <- cantau[g[, t]]

                canll <- 0.5 * sum(log(cantaug)) -
                         0.5 * quad.form(prec.cor, sqrt(cantaug) * res.t)

                R    <- canll - curll[t] +
                        dgamma(cantau[k], tau.alpha, tau.beta, log=TRUE) -
                        dgamma(tau[k, t], tau.alpha, tau.beta, log=TRUE) +
                        tryCatch(
                          {dgamma(tau[k, t], aaa, bbb, log=TRUE)}, 
                          warning = function(w) {
                            print(paste("knot", k, ", day", t))
                            print(paste("aaa =", aaa))
                            print(paste("bbb =", bbb))
                            print(paste("tau[k, t] =", tau[k, t]))
                            print(paste("g[, t] =", g[, t]))
                        }) -
                        dgamma(cantau[k], aaa, bbb, log=TRUE)
                
                if (!is.na(R)) { if (log(runif(1)) < R) {
                  acc.tau.ns[mh.idx] <- acc.tau.ns[mh.idx] + 1
                  acc.tau[k, t] <- acc.tau[k, t] + 1
                  tau[, t]  <- cantau
                  taug[, t] <- cantaug
                  curll[t]  <- canll
                }}
              }#end ifelse
              
          }#end k
        }  # end t
        
        }  # fi knots == 1
     
      }  # fi temporal
      
      # update tau.alpha and tau.beta
      a.star <- tau.beta.a + nt * nknots * tau.alpha
      b.star <- tau.beta.b + sum(tau)
      tau.beta <- rgamma(1, a.star, b.star)

      lll <- mmm <- seq(0.5, 10, 0.1)
      for (l in 1:length(lll)) {
        lll[l] <- sum(dgamma(tau, mmm[l], tau.beta, log=T))
      }
      tau.alpha <- sample(mmm, 1, prob=exp(lll - max(lll)))
      
      # update candidate tau
      for (t in 1:nt) { for (k in 1:nknots) {
        if ((att.tau[k, t] > 50) & (iter < (burn / 2))) {
          if (acc.tau[k, t] / att.tau[k, t] < 0.25) { mh.tau[k, t] <- mh.tau[k, t] * 0.8 }
          if (acc.tau[k, t] / att.tau[k, t] > 0.50) { mh.tau[k, t] <- mh.tau[k, t] * 1.2 }
          acc.tau[k, t] <- att.tau[k, t] <- 0
        }
      } }
      
      for (i in 1:length(mh.tau.ns)) {
        if ((att.tau.ns[i] > 50) & (iter < (burn / 2))) {
          if (acc.tau.ns[i] / att.tau.ns[i] < 0.25) { mh.tau.ns[i] <- mh.tau.ns[i] * 0.8 }
          if (acc.tau.ns[i] / att.tau.ns[i] > 0.50) { mh.tau.ns[i] <- mh.tau.ns[i] * 1.2 }
          acc.tau.ns[i] <- att.tau.ns[i] <- 0
        }
      }
    }  # fi method == t
    
    # update rho and nu and alpha
    if (rho.prior == "disc") {  # only update rho and alpha
      mu <- x.beta + z.alpha * zg
      res <- y - mu
      lll.rho <- rep(NA, length(rhos))
      
      # storage for possible covariance parts
      can.prec.cor <- array(NA, dim=c(ns, ns, length(rhos)))
      can.logdet.prec <- rep(NA, length(rhos))
      can.rss <- matrix(NA, nrow=length(rhos), ncol=nt)
      for (l in 1:length(rhos)) {
      	D <- 1 - alpha + alpha * C.values[, l]  # eigenvalues with alpha
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
            
      # update alpha
      att.alpha <- att.alpha + 1
    
      norm.alpha <- qnorm(alpha)
      can.norm.alpha <- rnorm(1, norm.alpha, mh.alpha)
      can.alpha <- pnorm(can.norm.alpha)
    
      can.C <- CorFx(d=d, alpha=can.alpha, rho=rho, nu=nu)
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
            dnorm(can.norm.alpha, mean=alpha.m, sd=alpha.s, log=T) - 
            dnorm(norm.alpha, mean=alpha.m, sd=alpha.s, log=T)
    
      if (!is.na(R)) { if (log(runif(1)) < R) {
        alpha <- can.alpha
        C <- can.C
        sd.mtx <- can.sd.mtx
        prec.cor <- can.prec.cor
        logdet.prec <- can.logdet.prec
        cur.rss <- can.rss
        acc.alpha <- acc.alpha + 1
      }}
    
      if ((att.alpha > 50) & (iter < (burn / 2))) {
        if (acc.alpha / att.alpha < 0.25) { mh.alpha <- mh.alpha * 0.8 }
        if (acc.alpha / att.alpha > 0.50) { mh.alpha <- mh.alpha * 1.2 }
        acc.alpha <- att.alpha <- 0
      }

    } else {
      att.rho <- att.rho + 1
      att.nu  <- att.nu + 1
      att.alpha <- att.alpha + 1
   
      logrho <- log(rho)
      can.logrho <- rnorm(1, logrho, mh.rho)
      can.rho <- exp(can.logrho)

      lognu  <- log(nu)
      if (!fixnu) {
        can.lognu <- rnorm(1, lognu, mh.nu)
      } else {
        can.lognu <- lognu
      }
      can.nu <- exp(can.lognu)
    
      norm.alpha <- qnorm(alpha)
      can.norm.alpha <- rnorm(1, norm.alpha, mh.alpha)
      can.alpha <- pnorm(can.norm.alpha)
    
      can.C <- CorFx(d=d, alpha=can.alpha, rho=can.rho, nu=can.nu)
      can.CC <- tryCatch(chol.inv(can.C, inv=T, logdet=T),
                         error = function(e) {
                           eig.inv(can.C, inv=T, logdet=T, mtx.sqrt=T)
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
            dnorm(lognu, lognu.m, lognu.s, log=T) +
            dnorm(can.norm.alpha, mean=alpha.m, sd=alpha.s, log=T) - 
            dnorm(norm.alpha, mean=alpha.m, sd=alpha.s, log=T)
    
      if (can.nu <= 10) {  # sometimes nu gets lost in the MCMC and ends up way too big
        if (!is.na(R)) { if (log(runif(1)) < R) {
          rho <- can.rho
          nu <- can.nu
          alpha <- can.alpha
          C <- can.C
          sd.mtx <- can.sd.mtx
          prec.cor <- can.prec.cor
          logdet.prec <- can.logdet.prec
          cur.rss <- can.rss
          acc.rho <- acc.rho + 1
          acc.nu  <- acc.nu + 1
          acc.alpha <- acc.alpha + 1
        }}
      }
    
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
    
      if ((att.alpha > 50) & (iter < (burn / 2))) {
        if (acc.alpha / att.alpha < 0.25) { mh.alpha <- mh.alpha * 0.8 }
        if (acc.alpha / att.alpha > 0.50) { mh.alpha <- mh.alpha * 1.2 }
        acc.alpha <- att.alpha <- 0
      }

    }
    
    if (skew) {
      # Skewness parameter
      vvv <- 1 / z.alpha.s^2
      mmm <- z.alpha.m
    
      for (t in 1:nt) {
      	taug.t <- sqrt(taug[, t])
      	prec.t <- sweep(prec.cor, 2, taug.t, "*") * taug.t 
        # prec.cov <- quad.form(prec.cor, diag(sqrt(taug[, t])))
        ttt <- zg[, t] %*% prec.t
        vvv <- vvv + ttt %*% zg[, t]
        mmm <- mmm + ttt %*% (y[, t] - x.beta[, t])
      }
    
      vvv <- 1 / vvv
      mmm <- vvv * mmm 
    
      z.alpha <- rnorm(1, mmm, sqrt(vvv))
    
      # z random effect
      mu <- y - x.beta - z.alpha * zg
      
      if (temporalz) {  # need to use MH sampling if there is a time series on the z terms
      	ts.z.update <- ts.sample.z(z.star=z.star, acc.z=acc.z, att.z=att.z, mh.z=mh.z, zg=zg,
      	                           phi=phi.z, att.phi=att.phi.z, acc.phi=acc.phi.z, mh.phi=mh.phi.z,
      	                           y=y, z.alpha=z.alpha, x.beta=x.beta, tau=tau, taug=taug, 
      	                           g=g, prec.cor=prec.cor)
      	
      	z.star <- ts.z.update$z.star
      	z     <- abs(z.star)
      	zg    <- ts.z.update$zg
      	phi.z <- ts.z.update$phi
      	att.z <- ts.z.update$att.z
      	acc.z <- ts.z.update$acc.z
      	att.phi.z <- ts.z.update$att.phi
      	acc.phi.z <- ts.z.update$acc.phi
      	
      	for (t in 1:nt) { 
          zg[, t] <- z[g[, t], t]
          for (k in 1:nknots) {
            if (att.z[k, t] > 50) {  # block accepting all knots for a day
              if (acc.z[k, t] / att.z[k, t] < 0.25) { mh.z[k, t] <- mh.z[k, t] * 0.8 }
              if (acc.z[k, t] / att.z[k, t] > 0.50) { mh.z[k, t] <- mh.z[k, t] * 1.2 }
              acc.z[k, t] <- att.z[k, t] <- 0
            }
          }
        } 
        
        if ((att.phi.z > 50) & (iter < (burn / 2)) ) {
          if (acc.phi.z / att.phi.z < 0.25) { mh.phi.z <- mh.phi.z * 0.8 }
          if (acc.phi.z / att.phi.z > 0.50) { mh.phi.z <- mh.phi.z * 1.2 }
          acc.phi.z <- att.phi.z <- 0
        }
      } else { 
        for (t in 1:nt) {
          taug.t <- sqrt(taug[,t])
          prec.t <- sweep(prec.cor, 2, taug.t, "*") * taug.t
          # prec.cov <- quad.form(prec.cor, diag(sqrt(taug[, t])))
          for (k in 1:nknots) {
            these <- which(g[, t] == k)
            r.1 <- y[these, t, drop=F] - x.beta[these, t, drop=F]
            r.2 <- y[-these, t, drop=F] - mu[-these, t, drop=F]
            prec.11 <- prec.t[these, these, drop=F]  # with cov
            prec.21 <- prec.t[-these, these, drop=F]
            lambda.l <- z.alpha^2 * sum(prec.11) + tau[k, t]
          
            mu.l <- z.alpha * sum(t(r.1) %*% prec.11 + t(r.2) %*% prec.21)
          
            e.z <- mu.l / lambda.l
            sd.z <- 1 / sqrt(lambda.l)
            z.new <- rTNorm(mn=e.z, sd=sd.z, lower=0, upper=Inf)
            if (z.new == Inf) {  # if z.new comes back Inf, then P(z > 0) = 0
              z.new = 0.00001
            }
        
            z[k, t] <- z.new
            zg[these, t] <- z.new
          }
        }
      } # fi temporal
    }    
  }  # end nthin
  
  # predictions
  if (predictions) {
  	if (cov.model == "matern") {
      s.11 <- alpha * matrix(matern(d11, rho, nu), np, np)
      s.12 <- alpha * matrix(matern(d12, rho, nu), np, ns)
    } else {
      s.11 <- alpha * matrix(exp(-d11 / rho), np, np)
      s.12 <- alpha * matrix(exp(-d12 / rho), np, ns)
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
      mup <- xp.beta + z.alpha * zgp + siggp * s.12.22.inv %*% (taug.t * res[, t])
      
      yp[, t] <- mup + siggp * t(corp.sd.mtx) %*% rnorm(np, 0, 1)
    }
   
  }
    
  # print(iter)
  # print("beta, tau.alpha, tau.beta")
  # print(c(beta, tau.alpha, tau.beta))
  # print("rho, alpha, nu")
  # print(c(rho, alpha, nu))

  # par(mfrow=c(3,2))
  # plot(y[,id1],pch=19,col=g[,id1],main=iter,ylim=range(y))
  # plot(y[,id2],pch=19,col=g[,id2],main=iter,ylim=range(y))
  # plot(1/sqrt(taug[,id1]),pch=19,col=g[,id1])
  # plot(1/sqrt(taug[,id2]),pch=19,col=g[,id2])
  # plot(s,col=g[,id1],ylim=0:1)
  # points(knots[,,id1],pch=19,col=1:nknots)
  # plot(s,col=g[,id2],ylim=0:1)
  # points(knots[,,id2],pch=19,col=1:nknots)
  
  keepers.tau[iter, , ]   <- tau
  keepers.beta[iter, ]    <- beta
  keepers.tau.alpha[iter] <- tau.alpha
  keepers.tau.beta[iter]  <- tau.beta
  keepers.rho[iter]       <- rho
  keepers.nu[iter]        <- nu
  keepers.alpha[iter]     <- alpha
  if (keep.knots & (nknots > 1)) {
    keepers.knots[iter, , , ] <- knots
  }
  
  if (predictions) {
    y.pred[iter, , ] <- yp
  }
  if (skew) {
    keepers.z.alpha[iter] <- z.alpha
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
  	acc.rate.alpha <- round(acc.alpha / att.alpha, 3)
  	
  	if (iter < burn) {
  	  begin <- max(1, (iter - 2000))
  	} else {
  	  begin <- burn
  	}
  	# print(mh.tau.ns)
    if (iterplot) {
      if (skew) {
        par(mfrow=c(3, 5))
      } else {
        par(mfrow=c(3, 4))
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
      # plot(keepers.beta[begin:iter, 1], type="l")
      plot(keepers.tau.alpha[begin:iter], type="l")
      plot(keepers.tau.beta[begin:iter], type="l")
      
      title.rho <- paste("acc =", acc.rate.rho)
      plot(keepers.rho[begin:iter], type="l", main=title.rho)
      title.nu <- paste("acc =", acc.rate.nu)
      plot(keepers.nu[begin:iter], type="l", main=title.nu)
      
      title.tau.1 <- paste("acc = ", acc.rate.tau[1, 1])
      title.tau.2 <- paste("acc = ", acc.rate.tau[1, 10])
      title.tau.3 <- paste("acc = ", acc.rate.tau[1, 21])
      plot(keepers.tau[begin:iter, 1, 1], type="l", main=title.tau.1)
      plot(keepers.tau[begin:iter, 1, 10], type="l", main=title.tau.2)
      plot(keepers.tau[begin:iter, 1, 21], type="l", main=title.tau.3)
      
      title.alpha <- paste("acc =", acc.rate.alpha)
      plot(keepers.alpha[begin:iter], type="l", main=title.alpha)
      if (skew) {
        plot(keepers.z.alpha[begin:iter], type="l")
        
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
    }
    
    toc <- proc.time()
    cat("\t iter", iter, "\n")
    # cat("\t elapsed time", (toc - tic)[3], "\n")
    tic <- proc.time()
    # cat("\t nu = ", nu, "\n")
    # cat("\t alpha = ", alpha, "\n")
    # cat("\t rho = ", rho, "\n")
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
  keepers.z.alpha <- NULL
} else {
  keepers.z <- keepers.z[return.iters, , ]
  keepers.z.alpha <- keepers.z.alpha[return.iters]
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
                alpha=keepers.alpha[return.iters],
                yp=y.pred,
                z.alpha=keepers.z.alpha,
                z=keepers.z,
                knots=keepers.knots,
                avgparts=keepers.avgparts,
                phi.z=keepers.phi.z,  # ts
                phi.w=keepers.phi.w,  # ts
                phi.tau=keepers.phi.tau  # ts
                )

return(results)
}#end mcmc()

