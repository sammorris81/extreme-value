#########################################################################
# MCMC 
#
# TODO: Add in model description here
#########################################################################

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL, 
                 thresh=0, thresh.quant=T, nknots=1, tau.by.knots=T,         
                 iters=5000, burn=1000, update=100, thin=1, scale=T,
                 iterplot=F, plotname=NULL,
                 # initial values
                 beta.init=NULL, tau.init=1,
                 tau.alpha.init=0.1, tau.beta.init=0.1,
                 rho.init=0.5, nu.init=0.5, alpha.init=0.5,
                 delta.init=0,
                 # priors
                 beta.m=0, beta.s=10, 
                 tau.alpha.m=0, tau.alpha.s=1, 
                 tau.beta.a=0.1, tau.beta.b=0.1,
                 logrho.m=-2, logrho.s=1,
                 lognu.m=-1.2, lognu.s=1,
                 alpha.m=0, alpha.s=1, # z.s=0.1,
                 # debugging settings
                 debug=F, knots.init, z.init, 
                 fixknots=F, fixz=F, fixbeta=F, 
                 fixtau=F, fixtau.alpha=F, fixtau.beta=F,
                 fixrho=F, fixnu=F, fixalpha=F,
                 fixdelta=F){
    
  start.time <- proc.time()
  
  ##############################################
  # Initial setup
  ##############################################
  ns <- nrow(y)  # number of sites
  nt <- ncol(y)  # number of days
    
  # rescale the x and y coordinates to be in [0, 1] x [0, 1]
  predictions <- !is.null(s.pred) & !is.null(x.pred)
  np <- 0
  y.pred <- NULL
  if (predictions) {
    np <- nrow(s.pred)
    if (scale) {
      s.unscale      <- s
      s.pred.unscale <- s.pred
      s.scale        <- ScaleLocs(rbind(s.unscale, s.pred.unscale))
      s              <- s.scale[1:ns, ]
      s.pred         <- s.scale[(ns + 1):(ns + np), ]
    }
    d12 <- rdist(s.pred, s)
    d11 <- rdist(s.pred, s.pred)
    diag(d11) <- 0
    y.pred <- array(0, c(np, nt, iters))
  } else {
    if (scale) {
      s.unscale <- s
      s <- ScaleLocs(s = s)
    }
  }    

  d       <- rdist(s)  # distance between sites
  diag(d) <- 0
    
  p <- dim(x)[3]  # number of covariates
  if (is.null(nknots)) { nknots = 1 }  # even if not using partition, still need
                                       # nknots
    
  # store the loc/day for observed values below thresh.
  if (thresh.quant & thresh > 0) {  # threshold based on sample quantiles
    thresh.data <- quantile(y, thresh, na.rm=T)
  } else {                        # threshold based on fixed value
    thresh.data <- thresh
  }
  thresh.mtx   <- matrix(thresh.data, ns, nt)  # need thresh in matrix form
  thresh.obs   <- !is.na(y) & (y <= thresh.mtx)
    
  # if there are some missing days, store the indices before
  # the first imputation. then set initial value as the average ozone
  # for all sites on the day.
  missing.obs <- is.na(y)
  for (t in 1:nt) {
    missing.sites       <- missing.obs[, t]
    y[missing.sites, t] <- mean(y[, t], na.rm=T)
  }
    
  # initialize partition
  knots <- vector(mode="list", length=nt)
  knots.min.1 <- range( s[, 1] )[1]
  knots.max.1 <- range( s[, 1] )[2]
  knots.min.2 <- range( s[, 2] )[1]
  knots.max.2 <- range( s[, 2] )[2]
  for (t in 1:nt) {
    knots[[t]]     <- matrix(NA, nknots, 2)
    knots[[t]][, 1] <- runif(n=nknots, min=knots.min.1, max=knots.max.1)
    knots[[t]][, 2] <- runif(n=nknots, min=knots.min.2, max=knots.max.2)
  }
    
  if(fixknots){knots <- knots.init}  # debug
       
  partition <- Membership(s=s, knots=knots)
  
  # initialize random effects
  z.knots <- matrix(1, nrow=nknots, ncol=nt)
  if (fixz) {z.knots <- z.init}
  z.sites <- ZBySites(z.knots, partition, nknots)
    
  # initialize parameters
  if (is.null(beta.init)) {
    beta <- rep(0, p)
    beta[1] <- mean(y)  # set initial intercept
  } else {
  	beta <- beta.init
  }
  if(fixbeta){beta <- beta.init}  #debug
   
  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
  	x.beta[, t] <- x[, t, ] %*% beta
  }
  
  if (length(tau.init) == 1) {
  	cat("  initializing all tau terms as", tau.init, "\n")
    tau.knots <- matrix(data=tau.init, nrow=nknots, ncol=nt)
  } else {
    if (nrow(tau.init) != nknots | ncol(tau.init) != nt) {
      stop("tau.init must be a matrix with nknots rows and nt cols")
    }
    tau.knots <- matrix(data=tau.init, nrow=nknots, ncol=nt)
  }
  
  tau.sites <- TauSites(tau.knots=tau.knots, partition=partition, 
                            nknots=nknots)
  
  # easier to keep calculations in the precision scale for MCMC
  sigma.knots <- 1 / tau.knots
  sigma.sites <- 1 / tau.sites
  tau.alpha <- tau.alpha.init
  tau.beta  <- tau.beta.init
  
  # initialize spatial covariance
  rho <- rho.init
  logrho <- log(rho)
  nu  <- nu.init
  lognu <- log(nu)
  alpha <- alpha.init
    
  cor.mtx     <- SpatCor(d=d, alpha=alpha, rho=rho, nu=nu)
  prec.cor    <- cor.mtx$prec.cor
  logdet.prec <- cor.mtx$logdet.prec
  
  if (delta.init < -1 | delta.init > 1) {
    stop("delta.init must be between -1 and 1.")
  } else {
    delta <- delta.init
  }
  
  # MH tuning params
  acc.w <- att.w <- mh.w <- rep(0.1, nt)  # knot locations
  acc.tau <- att.tau <- mh.tau <- matrix(1, nrow=nknots, ncol=nt)
  # acc.sigma.alpha <- att.sigma.alpha <- mh.sigma.alpha <- 0.1 
  acc.delta <- att.delta <- mh.delta <- 0.01  
  acc.rho   <- att.rho   <- mh.rho   <- 1
  acc.nu    <- att.nu    <- mh.nu    <- 1
  acc.alpha <- att.alpha <- mh.alpha <- 1
  
  # storage
  if (debug) { print("storage") }
  keepers.z <- keepers.tau <- array(NA, dim=c(iters, nknots, nt))
  keepers.beta <- matrix(NA, nrow=iters, ncol=p)
  keepers.ll <- matrix(NA, nrow=iters, ncol=nt)
  keepers.tau.alpha <- keepers.tau.beta <- rep(NA, iters)
  keepers.delta <- keepers.rho <- keepers.nu <- keepers.alpha <- rep(NA, iters)
  
  # initial values
  mu <- x.beta + delta * z.sites
  cur.ll <- LLike(y=y, x.beta=x.beta, tau.sites=tau.sites, delta=delta, prec.cor=prec.cor, 
                  logdet.prec=logdet.prec, z.sites=z.sites, log=T)
  
  for (iter in 1:iters) { for (ttt in 1:thin) {
    
    # impute data below threshold
    if (debug) { print("impute") }
    if (thresh != 0) {
      mu <- x.beta + delta * z.sites
      thresh.mtx.fudge <- 0.99999 * thresh.mtx  # numerical stability
      y.imputed <- matrix(y, ns, nt)
      
      cor <- matern(u=d, phi=rho, kappa=nu)  # doesn't change each day
      for (t in 1:nt) {
      	these.thresh.obs <- thresh.obs[, t]
      	these.missing.obs <- missing.obs[, t]
        
        # spatial error
        theta.t <- sqrt((1 - delta^2) * alpha / tau.sites[, t]) * t(chol(cor)) %*% 
                   rnorm(ns, 0, 1)
        
        # nugget error
        # new expected value and standard deviation
        e.y <- mu[, t] + theta.t
        s.y <- sqrt((1 - delta^2) * (1 - alpha) / tau.sites[, t])
        upper.y <- thresh.mtx[, t]
        
        y.impute.t <- rTNorm(mn=e.y, sd=s.y, lower=-Inf, upper=upper.y)
        
        # if any y.impute.t come back as -Inf it's because P(Y < T) = 0
        usefudge <- y.impute.t == -Inf
        y.impute.t[usefudge] <- thresh.mtx.fudge[usefudge, t]
        # y.impute.t <- ifelse(y.impute.t == -Inf, thresh.mtx.fudge[, t], y.impute.t)
        y.imputed[these.thresh.obs, t] <- y.impute.t[these.thresh.obs]
        
        y.missing.t <- rnorm(n=ns, mean=e.y, sd=s.y)
        y.imputed[these.missing.obs, t] <- y.missing.t[these.missing.obs]
      }  # end nt
      
      y <- y.imputed
      
    }  # fi thresh != 0
    
    # update random effects
    if (debug) {print("random z")}
    if (!fixz) {  # debug
      for (t in 1:nt) {
        if (nknots > 1) {
          for (k in 1:nknots) {
            these    <- which(partition[, t] == k)
            n.these  <- length(these)
            r1       <- y[these, t] - x.beta[these, t]
            r2       <- y[-these, t] - x.beta[-these, t] - delta * z.sites[-these, t]
            
            # vvv      <- sigma[, t] * (1 - delta^2)
            prec.t     <- diag(tau.sites[, t]) %*% prec.cor %*% diag(tau.sites[, t])
            prec.11.t  <- prec.t[these, these]
            prec.21.t  <- prec.t[-these, these]
            
            mu.z     <- delta * sum(t(r1) %*% prec.11.t + t(r2) %*% prec.21.t)
            lambda.z <- delta^2 * sum(prec.11.t) + (1 - delta^2) * tau.knots[k, t]
            mn.z     <- mu.z / lambda.z
            sd.z     <- sqrt((1 - delta^2) / lambda.z)
        
            z.new   <- rTNorm(mn=mn.z, sd=sd.z, lower=0, upper=Inf) 
            # if any z come back as Inf it's because P(Y > T) = 0
            if (z.new == Inf) {
              z.new = 0.0001
            }

            z.knots[k, t]     <- z.new
            z.sites[these, t] <- z.new
          }  # end nknots
        } else if (nknots == 1) {
          prec.cov <- diag(sqrt(tau.sites[, t])) %*% prec.cor %*% diag(sqrt(tau.sites[, t]))
          r        <- y[, t] - x.beta[, t]
          mu.z     <- delta * sum(prec.cov %*% r)
          lambda.z <- delta^2 * sum(prec.cov) + (1 - delta^2) * tau.knots[, t]
          mn.z     <- mu.z / lambda.z
          sd.z     <- sqrt((1 - delta^2) / lambda.z)
          z.new    <- rTNorm(mn=mn.z, sd=sd.z, lower=0, upper=Inf)
          
          # if any z come back as Inf it's because P(Y > T) = 0
          if (z.new == Inf) {
            z.new = 0.0001
          }
          
          z.knots[, t] <- z.new
          z.sites[, t] <- z.new
        }  # fi nknots
      }  # end nt
    }  # fi !fixz
    
    # update means
    mu <- x.beta + delta * z.sites
    
    # update partitions
    if (debug) { print("knots") }
    if (!fixknots) {  # debug
      if (nknots > 1) {
        can.knots <- vector(mode="list", length=1)
        for (t in 1:nt) {
          can.knots[[1]] <- matrix(NA, nrow=nknots, ncol=2)
          att.w[t] <- att.w[t] + 1
          res.t <- matrix((y[, t] - mu[, t]), ns, 1)
          tau.sites.t <- matrix(tau.sites[, t], ns, 1)
          ss.t <- SumSquares(res.t, prec.cor, tau.sites.t) / (1 - delta^2)
          cur.ll <- 0.5 * sum(log(tau.sites.t)) - 0.5 * ss.t
          
          for (k in 1:nknots) {
        	can.knots[[1]][k, ] <- rnorm(2, knots[[t]][k, ], mh.w[t])
          }
        
          can.partition.t <- Membership(s=s, knots=can.knots)
          z.knots.t <- matrix(z.knots[, t], nknots, 1)
          can.z.sites.t <- ZBySites(z.knots.t, can.partition.t, nknots)
          can.tau.sites.t <- TauSites(tau.knots, can.partition.t, nknots)
          can.mu.t <- x.beta[, t] + delta * can.z.sites.t
          can.res.t <- matrix((y[, t] - can.mu.t), ns, 1)
          can.ss.t <- SumSquares(can.res.t, prec.cor, can.tau.sites.t) / (1 - delta^2)
          
          if (length(can.tau.sites.t) != ns) {
            stop("can.tau.sites.t is the wrong length.")
          }
          
          if (length(can.ss.t) != 1) {
            stop("can.ss.t is the wrong length")
          }

          can.ll <- 0.5 * sum(log(can.tau.sites.t)) - 0.5 * can.ss.t
          
          rej <- sum(can.ll - cur.ll)  # prior is uniform and candidate is symmetric
          
          if (!is.na(rej)) {if (-rej < rexp(1, 1)) {
            knots[[t]]     <- can.knots[[1]]
            partition[, t] <- can.partition.t
            z.sites[, t]   <- can.z.sites.t
            tau.sites[, t] <- can.tau.sites.t
            mu[, t]        <- can.mu.t
            acc.w[t]       <- acc.w[t] + 1
            cur.ll         <- can.ll
          }}
          
        } # end nt
        
        partition <- Membership(s=s, knots=knots)
        z.sites <- ZBySites(z.knots, partition, nknots=nknots)
        tau.sites <- TauSites(tau.knots, partition, nknots)
        mu <- x.beta + delta* z.sites
      }  # fi nknots > 1
    }  # fi !fixknots
    
    # update beta
    if (debug) { print("beta") }
    if (!fixbeta) {  # debug
      vvv <- diag(p) / beta.s^2
      mmm <- rep(beta.m, p)
      
      for (t in 1:nt) {
        x.t <- x[, t, ]
        prec.t <- sqrt(diag(tau.sites[, t])) %*% prec.cor %*% sqrt(diag(tau.sites[, t]))
        ttt <- t(x.t) %*% prec.t / (1 - delta^2)
        vvv <- vvv + ttt %*% x.t
        mmm <- mmm + ttt %*% (y[, t] - delta * z.sites[, t])
      }
      
      vvv <- solve(vvv)
      beta <- vvv %*% mmm + t(chol(vvv)) %*% rnorm(p)
      
      for (t in 1:nt) {
  	    x.beta[, t] <- x[, t, ] %*% beta
      }
    }  #fi !fixbeta
    
    mu  <- x.beta + delta * z.sites
    
    if (debug) { print("delta") }
    if (!fixdelta) {
      att.delta <- att.delta + 1
      
      res <- y - mu
      cur.rss <- SumSquares(res, prec.cor, tau.sites)
      # delta.tan <- tan(delta * pi / 2)
      
      can.delta <- rnorm(1, delta, mh.delta)
      # can.delta.tan <- rnorm(1, delta.tan, mh.delta)
      
      if (can.delta > -1 & can.delta < 1) {
      	# can.delta <- atan(can.delta.tan) * 2 / pi
	    can.res <- y - (x.beta + can.delta * z.sites)
	    can.rss <- SumSquares(can.res, prec.cor, tau.sites)

	    rej <- -0.5 * nt * ns * (log(1 - can.delta^2) - log(1 - delta^2)) - 
	           0.5 * (sum(can.rss) / (1 - can.delta^2) - sum(cur.rss) / (1 - delta^2))
	    
	    if (length(rej) > 1) {
	      stop("rej for delta is too long")
	    }
	           
	    if (!is.na(rej)) { if (-rej < rexp(1, 1)) {
	      delta     <- can.delta
	      acc.delta <- acc.delta + 1
	    }}
	  }  # fi can.delta
      
      # lll <- mmm <- seq(-0.9999, 0.9999, 0.01)
      # for (l in 1:length(lll)) {
      	# lll[l] <- sum(LLike(y, x.beta, tau.sites, mmm[l], prec.cor, logdet.prec, z.sites))
      # }
      # delta <- sample(mmm, 1, prob=exp(lll - max(lll)))
      
    }  # fi !fixdelta
    
    #### Spatial correlation
    
    mu  <- x.beta + delta * z.sites
    res <- y - mu
    
    # update tau (inverse-sill)
    if (debug) { print("tau") }
    if (!fixtau) {  # tau.by.knots=F means same tau for all knots
      if (!tau.by.knots | (nknots == 1)) {
      # if (!tau.by.knots) {
      	rss <- rep(NA, nt)
      	for (t in 1:nt) {
      	  rss[t] <- t(res[, t]) %*% prec.cor %*% res[, t] / (1 - delta^2)
      	}
        alpha.star  <- tau.alpha + nknots / 2 + ns / 2
        beta.star   <- tau.beta + colSums(z.knots^2) / 2 + rss / 2
        tau.knots <- rgamma(nt, alpha.star, beta.star)
        tau.knots <- matrix(tau.knots, nknots, nt)
      } else {
        for (k in 1:nknots) {
          att.tau[k, ] <- att.tau[k, ] + 1
          can.tau.knots <- tau.knots
          can.logtau.knots <- logtau.knots <- log(tau.knots)
          can.logtau.knots[k, ] <- rnorm(nt, logtau.knots[k, ], mh.tau[k, ])
          can.tau.knots[k, ] <- exp(can.logtau.knots[k, ])
          can.tau.sites <- TauSites(can.tau.knots, partition, nknots)
                    
          can.rss <- SumSquares(res, prec.cor, can.tau.sites)
          cur.rss <- SumSquares(res, prec.cor, tau.sites)
  
          ns.k <- colSums(partition == k)
  
          if ((nrow(can.tau.sites) != ns) | (ncol(can.tau.sites) != nt)) {
            stop("can.tau.sites is not the correct dimensions")
          }

          if ((length(can.rss) != nt) | (length(cur.rss) != nt) | (length(ns.k) != nt)) {
            stop("one of can.rss, cur.rss, or ns.k is not the correct length")
          }  
  
          # can.ll <- LLike(y, x.beta, can.tau.sites, delta, prec.cor, logdet.prec, z.sites, log=T)
          # cur.ll <- LLike(y, x.beta, tau.sites, delta, prec.cor, logdet.prec, z.sites, log=T)
          
          # can.ll <- 0.5 * (ns.k + 1) * can.logtau.knots[k, ] - 
                    # 0.5 * (z.knots[k, ]^2 * can.tau.knots[k, ] + can.rss / (1 - delta^2))
            
          # cur.ll <- 0.5 * (ns.k + 1) * logtau.knots[k, ] - 
                    # 0.5 * (z.knots[k, ]^2 * tau.knots[k, ] + cur.rss / (1 - delta^2))
          
          tau.alpha.star <- tau.alpha + 0.5 + 0.5 * ns.k
          tau.beta.star  <- tau.beta + 0.5 * z.knots[k, ]^2
          
          rej <- dgamma(can.tau.knots[k, ], tau.alpha.star, tau.beta.star, log=T) -
                 dgamma(tau.knots[k, ], tau.alpha.star, tau.beta.star, log=T) - 
                 0.5 * (can.rss - cur.rss) / (1 - delta^2)
                 
  
          if (length(rej) != nt) {
            stop("rej is not the correct length")
          }

          # accept / reject       
          # mh.compare <- rexp(nt, 1)
          # mh.update <- (-rej < mh.compare) & (!is.na(rej))
          # tau.knots[k, mh.update] <- can.tau.knots[k, mh.update]
          # acc.tau[k, mh.update] <- acc.tau[k, mh.update] + 1
          for (t in 1:nt) {
            if (!is.na(rej[t])) { if (-rej[t] < rexp(1, 1)) {
              tau.knots[k, t] <- can.tau.knots[k, t]
              acc.tau[k, t] <- acc.tau[k, t] + 1
            }}
          }
          
          tau.sites <- TauSites(tau.knots, partition, nknots)
        }  # end nknots
      }  # fi nknots == 1
      # cat("tau.knots", tau.knots, "\n")
      
    }  # fi !fixtau
    
    # update tau.alpha and tau.beta
    if (debug) { print("tau.alpha and tau.beta") }
        
    if (!fixtau.beta) {
      a.star <- tau.beta.a + nt * nknots * tau.alpha
      b.star <- tau.beta.b + sum(tau.knots)
      
      tau.beta <- rgamma(1, a.star, b.star)
    }

#     Update the sigma.alpha terms on a grid
#     lll<-mmm<-seq(.5,30,.5)
#     for(l in 1:length(lll)){
#       lll[l]<-sum(dgamma(1/r,mmm[l],sig.r,log=T))
#     }
#     xi.r<-sample(mmm,1,prob=exp(lll-max(lll)))
    
    if (!fixtau.alpha) {
      lll <- mmm <- seq(0.5, 30, 0.5)
      for (l in 1:length(lll)) {
      	ll.tau <- dgamma(tau.knots, mmm[l], tau.beta, log=T)
        lll[l] <- sum(ll.tau)
      }
      # cat("lll is", lll, "\n")
      tau.alpha <- sample(mmm, 1, prob=exp(lll - max(lll)))
      
    }  # fi !fixtau.alpha
              
    # rho and nu
    if (debug) { print("rho and nu") }
    if (!fixrho | !fixnu) {
      
      cur.rss <- SumSquares(res, prec.cor, tau.sites)
      
      if (!fixrho) {
      	att.rho    <- att.rho + 1
        logrho     <- log(rho)
        can.logrho <- rnorm(1, logrho, mh.rho)
      } else {
        can.logrho <- logrho <- log(rho)
      }
      can.rho <- exp(can.logrho)
      
      if (!fixnu) {
      	att.nu    <- att.nu + 1
        lognu     <- log(nu)
        can.lognu <- rnorm(1, lognu, mh.nu)
      } else {
        can.lognu <- lognu <- log(nu)
      }
      can.nu <- exp(can.lognu)
      
      # if (!fixalpha) {
        # att.alpha     <- att.alpha + 1
        # phi.alpha     <- qnorm(alpha)
        # can.phi.alpha <- rnorm(1, phi.alpha, mh.alpha)
      # } else {
        # can.phi.alpha <- phi.alpha <- qnorm(alpha) 
      # }
      # can.alpha <- pnorm(can.phi.alpha)
      
      if (can.nu <= 10 & can.rho < 3) {  # for numerical stability
        # can.cor.mtx <- SpatCor(d=d, alpha=can.alpha, rho=can.rho, nu=can.nu)
        can.cor.mtx     <- SpatCor(d=d, alpha=alpha, rho=can.rho, nu=can.nu)
        can.prec.cor    <- can.cor.mtx$prec.cor
        can.logdet.prec <- can.cor.mtx$logdet.prec
        
        can.rss     <- SumSquares(res, can.prec.cor, tau.sites)    
                      
        rej <- dnorm(can.logrho, logrho.m, logrho.s, log=T) - 
               dnorm(logrho, logrho.m, logrho.s, log=T) +
               dnorm(can.lognu, lognu.m, lognu.s, log=T) -
               dnorm(lognu, lognu.m, lognu.s, log=T) - 
               0.5 * sum(can.rss - cur.rss) / (1 - delta^2) + 
               0.5 * nt * (can.logdet.prec - logdet.prec) # + 
               # dnorm(can.phi.alpha, mean=alpha.m, sd=alpha.s, log=T) - 
               # dnorm(phi.alpha, mean=alpha.m, sd=alpha.s, log=T)
        
        if (!is.na(rej)) { if (-rej < rexp(1, 1)) {
          rho         <- can.rho
          nu          <- can.nu
          # alpha       <- can.alpha
          prec.cor    <- can.prec.cor
          logdet.prec <- can.logdet.prec
          if (!fixrho) { acc.rho <- acc.rho + 1 }
          if (!fixnu) { acc.nu <- acc.nu + 1 }
          # if (!fixalpha) { acc.alpha <- acc.alpha + 1 }
        }} 
      } # fi can.nu <= 10
    } # fi !fixrho | !fixnu
    
    # alpha
    if (debug) { print("alpha") }
    if (!fixalpha) {
    
      cur.rss  <- SumSquares(res, prec.cor, tau.sites)
                      
      att.alpha      <- att.alpha + 1
      norm.alpha     <- qnorm(alpha)
      can.norm.alpha <- rnorm(1, norm.alpha, mh.alpha)
      can.alpha      <- pnorm(can.norm.alpha)
      
      can.cor.mtx     <- SpatCor(d=d, alpha=can.alpha, rho=rho, nu=nu)
      can.prec.cor    <- can.cor.mtx$prec.cor
      can.logdet.prec <- can.cor.mtx$logdet.prec
    
      can.rss <- SumSquares(res, can.prec.cor, tau.sites)
                      
      rej <- dnorm(can.norm.alpha, mean=alpha.m, sd=alpha.s, log=T) - 
             dnorm(norm.alpha, mean=alpha.m, sd=alpha.s, log=T) - 
             0.5 * sum(can.rss - cur.rss) / (1 - delta^2) + 
             0.5 * nt * (can.logdet.prec - logdet.prec)

      if (!is.na(rej)) { if (-rej < rexp(1, 1)) {
        alpha       <- can.alpha
        prec.cor    <- can.prec.cor
        logdet.prec <- can.logdet.prec
        acc.alpha   <- acc.alpha + 1
      }}

    }  # fi !fixalpha
    
  }  # end nthin
  
  ##############################################
  # Spatial Predictions
  ##############################################

  if (np > 0) { 	
  	mu <- x.beta + delta * z.sites  # just want to make sure it's up to date
    partition.pred <- Membership(s=s.pred, knots=knots)
    z.sites.pred   <- ZBySites(z.knots=z.knots, partition=partition.pred, nknots)
    tau.sites.pred <- TauSites(tau.knots=tau.knots, partition=partition.pred, nknots)
    res.std <- (y - mu) * sqrt(tau.sites)
    
    # these don't change from day to day
    s.11 <- diag(1, nrow=np)
    s.12 <- matrix(alpha * matern(d12, phi=rho, kappa=nu), nrow=np, ncol=ns)
    s.12.22.inv <- s.12 %*% prec.cor
    
    for (t in 1:nt) {
      x.beta.pred.t    <- x.pred[, t, ] %*% beta
      z.sites.pred.t   <- z.sites.pred[, t]
      tau.sites.pred.t <- tau.sites.pred[, t]
      tau.sites.t      <- tau.sites[, t]
      res.std.t        <- res.std[, t]
     
      sigma.sites.pred.t <- 1 / tau.sites.pred.t
      sigma.sites.t      <- 1 / tau.sites.t
      
      mu.pred.t <- x.beta.pred.t + delta * z.sites.pred.t
      
      # double check here. how do the tau sites fit in?
      # would this also work as
      # s.11 <- diag(1, nrow=np)
      # s.12 <- matrix(alpha * matern(d12, phi=rho, kappa=nu), nrow=np, ncol=ns)
      # s.22.inv <- prec.cor
      
      # s.12.22.inv <- s.12 %*% prec.cor
      e.y.pred <- mu.pred.t - sqrt(sigma.sites.pred.t) * s.12.22.inv %*% res.std.t
      v.y.pred <- (1 - delta^2) * (s.11 - s.12.22.inv %*% t(s.12))
      sd.y.pred <- sqrt(diag(v.y.pred)) * sqrt(sigma.sites.pred.t)
      
      if (length(e.y.pred) != np) {
        stop("e.y.pred is the wrong length")
      } else if (length(sd.y.pred) != np) {
      	stop("sd.y.pred is the wrong length")
      }
      
      y.pred[, t, iter] <- rnorm(np, mean=e.y.pred, sd=sd.y.pred)
      
    }  # end t
  }  # fi np > 0

  cur.ll <- LLike(y=y, x.beta=x.beta, tau.sites=tau.sites, delta=delta, 
                  prec.cor=prec.cor, logdet.prec=logdet.prec, z.sites=z.sites, log=T)
  
  ##############################################
  # Keep track of iterations
  ##############################################
  if (debug) { print("keepers") }
  keepers.z[iter, , ]     <- z.knots
  keepers.beta[iter, ]    <- beta
  keepers.tau[iter, , ]   <- tau.knots
  keepers.tau.alpha[iter] <- tau.alpha
  keepers.tau.beta[iter]  <- tau.beta
  keepers.delta[iter]     <- delta
  keepers.rho[iter]       <- rho
  keepers.nu[iter]        <- nu
  keepers.alpha[iter]     <- alpha
  keepers.ll[iter, ]      <- cur.ll

  ##############################################
  # Display current value
  ##############################################
  if (debug) { print("plotting") }
  
  if (iterplot) {
    plotnow <- F
    # different behavior if running on server vs testing
    if (((iter %% update) == 0) && is.null(plotname)) {
  	  plotnow <- T
  	  plotmain <- ""
  	} else if ((iter == iters) && !is.null(plotname)) {
  	  plotnow <- T
  	  plotmain <- plotname
      plotname.file <- paste("plots/", plotname, sep="")
  	  pdf(file=plotname.file)
  	}
  	
  	if (plotnow) {  
  	  accrate.w     <- round(acc.w / att.w, 3)
  	  accrate.delta <- round(acc.delta / att.delta, 3)
  	  accrate.rho   <- round(acc.rho / att.rho, 3)
  	  accrate.nu    <- round(acc.nu / att.nu, 3)
  	  accrate.alpha <- round(acc.alpha / att.alpha, 3)
  	  accrate.tau   <- round(acc.tau / att.tau, 3)
  	  	
  	  par(mfrow=c(3, 4))
  	  # if (iter > burn) {
  	    # start <- burn
  	  # } else {
  	    # start <- 1
  	  # }
  	  start <- 1
  	  cat("dim(y.pred) is ", dim(y.pred))
  	  plot(keepers.beta[start:iter, 1], ylab="beta0", xlab="iteration", 
           type="l")
      # plot(keepers.beta[start:iter, 2], ylab="beta1", xlab="iteration",
           # type="l", main=plotmain)
      # plot(keepers.beta[start:iter, 3], ylab="beta2", xlab="iteration",
           # type="l")
      plot(keepers.ll[start:iter], ylab="loglike", xlab="iteration",
           type="l")
      # plot(keepers.z[start:iter, 1, 8], ylab="z 1, 8", xlab="iteration", 
           # type="l")
      # plot(keepers.z[start:iter, 3, 6], ylab="z 3, 6", xlab="iteration",
           # type="l")
      # plot(keepers.z[start:iter, 3, 16], ylab="z 3, 16", xlab="iteration", 
           # type="l")
      # plot(keepers.z[start:iter, 1, 20], ylab="z 1, 20", xlab="iteration",
           # type="l")
      plot(y.pred[1, 1, start:iter], ylab="ypred 1, 1", xlab="iteration", 
           type="l")
      plot(y.pred[5, 1, start:iter], ylab="ypred 5, 1", xlab="iteration", 
           type="l")
      plot(keepers.delta[start:iter], ylab="delta", xlab="iteration", 
           type="l", main=bquote("ACCR" == .(accrate.delta)))
      plot(keepers.rho[start:iter], ylab="rho", xlab="iteration", 
            type="l", main=bquote("ACCR" == .(accrate.rho)))
      plot(keepers.nu[start:iter], ylab="nu", xlab="iteration", 
           type="l", main=bquote("ACCR" == .(accrate.nu)))
      plot(keepers.alpha[start:iter], ylab="alpha", xlab="iteration", 
           type="l", main=bquote("ACCR" == .(accrate.alpha)))
      # plot(keepers.tau[start:iter, , 1], ylab="sigma 1", xlab="iteration", 
           # type="l")
      # plot(keepers.tau[start:iter, , 3], ylab="sigma 3", xlab="iteration",
           # type="l")
      # plot(keepers.tau[start:iter, 1, 1], ylab="tau 1, 1", xlab="iteration", 
           # type="l", main=bquote("ACCR" == .(accrate.tau[1, 1])))
      # plot(keepers.tau[start:iter, 1, 20], ylab="tau 1, 20", xlab="iteration", 
           # type="l", main=bquote("ACCR" == .(accrate.tau[1, 20])))
      plot(keepers.tau[start:iter, 1, 1], ylab="tau 1, 1", xlab="iteration", 
           type="l")
      plot(keepers.tau[start:iter, 1, 20], ylab="tau 1, 20", xlab="iteration", 
           type="l")
      plot(keepers.tau.alpha[start:iter], ylab="tau.alpha", xlab="iteration",
           type="l")
      plot(keepers.tau.beta[start:iter], ylab="tau.beta", xlab="iteration",
           type="l")
  	}
  	
  	if ((iter == iters) && !is.null(plotname)) {
  	  dev.off()
  	}
  }  # fi iterplot
  
  ##############################################
  # MH adapt candidate distributions
  ##############################################
  if (debug) { print("MH update") }
  if (iter < burn / 2) { 
  	
  	accrate.w <- acc.w / att.w
    for (t in 1:nt) {
      if (att.w[t] > 50) {
        if (accrate.w[t] < 0.25) { mh.w[t] <- mh.w[t] * 0.8 }
        if (accrate.w[t] > 0.50) { mh.w[t] <- mh.w[t] * 1.2 }
        acc.w[t] <- att.w[t] <- 0.1
      }
    }
    
    accrate.tau <- acc.tau / att.tau
    for (t in 1:nt) {
      for (k in 1:nknots) {
        if (att.tau[k, t] > 50) {
          if (accrate.tau[k, t] < 0.25) { 
            mh.tau[k, t] <- mh.tau[k, t] * 0.8 
          }
          if (accrate.tau[k, t] > 0.50) {
            mh.tau[k, t] <- mh.tau[k, t] * 1.2
          }
        } 
      }
    }
    #cat("accrate.tau is", accrate.tau)
    
    # if (att.sigma.alpha > 50) {
      # if (acc.sigma.alpha / att.sigma.alpha < 0.25) { 
        # mh.sigma.alpha <- mh.sigma.alpha * 0.8
      # }
      # if (acc.sigma.alpha / att.sigma.alpha > 0.50) {
        # mh.sigma.alpha <- mh.sigma.alpha * 1.2
      # }
    # }
    
    accrate.delta <- acc.delta / att.delta
    if (att.delta > 50) {
      cat("mh.delta is", mh.delta, "\n")
      cat("accrate.delta is", accrate.delta, "\n")
      if (accrate.delta < 0.25) { mh.delta <- mh.delta * 0.8 }
      if (accrate.delta > 0.50) { mh.delta <- mh.delta * 1.2 }
      acc.delta <- att.delta <- 0.1
    }
    
    accrate.rho <- acc.rho / att.rho
    if (att.rho > 50) {
      if (accrate.rho < 0.25) { mh.rho <- mh.rho * 0.8 }
      if (accrate.rho > 0.50) { mh.rho <- mh.rho * 1.2 }
      acc.rho <- att.rho <- 0.1
      # print(mh.rho)
    }
    
    accrate.nu <- acc.nu / att.nu
    if (att.nu > 50) {
      if (accrate.nu < 0.25) { mh.nu <- mh.nu * 0.8 }
      if (accrate.nu > 0.50) { mh.nu <- mh.nu * 1.2 }
      acc.nu <- att.nu <- 0.1
      # print(mh.nu)
    }
    
    accrate.alpha <- acc.alpha / att.alpha
    if (att.alpha > 50) {
      if (accrate.alpha < 0.25) { mh.alpha <- mh.alpha * 0.8 }
      if (accrate.alpha > 0.50) { mh.alpha <- mh.alpha * 1.2 }
      acc.alpha <- att.alpha <- 0.1
      # print(mh.alpha)
    }
  }  # fi iter < burn / 2
  
  if ((iter %% update) == 0) {
    cat("    Iter", iter, "complete. \n")
  }
  if (debug) { print(paste("iter: ", iter)) }

  }  # end iter   
  
  ##############################################
  # Return output
  ##############################################
  stop.time <- proc.time()
  results   <- list(time=stop.time-start.time,
                    z=keepers.z,
                    beta=keepers.beta,
                    tau=keepers.tau,
                    tau.alpha=keepers.tau.alpha,
                    tau.beta=keepers.tau.beta,
                    delta=keepers.delta,
                    rho=keepers.rho,
                    nu=keepers.nu,
                    alpha=keepers.alpha,
                    yp=y.pred)
  
  return(results)
}
