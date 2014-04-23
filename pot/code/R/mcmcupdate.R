#########################################################################
# MCMC 
#
# TODO: Add in model description here
#########################################################################

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL, 
                 thresh=0, thresh.quant=T, 
                 r.model = "gamma", #also allow "fixed"
                 nknots=1,          
                 iters=5000, burn=1000, update=100, thin=1, scale=T,
                 iterplot=F, plotname=NULL,
                 # initial values
                 knots.init, z.init, beta.init=NULL, sigma.init=1,
                 rho.init=0.5, nu.init=0.5, alpha.init=0.5,
                 delta.init=0,
                 # priors
                 beta.m=0, beta.s=10, sigma.a=0.1, sigma.b=10,
                 # debugging settings
                 debug=F, 
                 fixknots=F, fixz=F, fixbeta=F, fixsigma=F
                 fixrho=F, fixnu=F, fixalpha=F,
                 fixdelta=F){
    
  start.time <- proc.time()
   
  ##############################################
  # Initial setup
  ##############################################
  
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
  if (thresh.quant) {               # threshold based on sample quantiles
    thresh.data <- quantile(y, thresh)
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
  z.sites <- matrix(NA, nrow=ns, ncol=nt)
  for (t in 1:nt) {
  	z.sites[, t] <- ZBySites(z.knots, partition)
  }
  
  if (fixz) {z <- matrix(z.init, nrow=nknots, ncol=nt)}
  
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
  
  # initialize spatial covariance
  rho <- rho.init
  logrho <- log(rho)
  nu  <- nu.init
  lognu <- log(nu)
  alpha <- alpha.init
  
  if (length(sigma) == 1) {
    sigma <- rep(sigma.init, nt)
  } else {
    sigma <- sigma.init
  }
  
  cor.mtx <- SpatCor(d=d, alpha=alpha, rho=rho, nu=nu)
  sig     <- cor.mtx$sig
  prec    <- cor.mtx$prec
  log.det <- cor.mtx$log.det
  
  if (delta.init < -1 || delta.init > 1) {
    stop("delta.init must be between -1 and 1.")
  } else {
    delta <- delta.init
  }
  
  # MH tuning params
  acc.w <- att.w <- mh.w <- rep(0.1, nt)  # knot locations
  acc.delta <- att.delta <- mh.delta <- 0.1  
  acc.alpha <- att.alpha <- mh.alpha <- 0.1
  acc.rho   <- att.rho   <- mh.rho   <- 0.1
  acc.nu    <- att.nu    <- mh.nu    <- 0.1
  
  for (iter in 1:iters) { for (ttt in 1:nthin) {
    
    # impute data below threshold
    if (debug) {print("impute")}
    if (thresh != 0) {
      mu <- x.beta + delta * z.sites
      thresh.mtx.fudge <- 0.99999 * thresh.mtx.k  # numerical stability
      y.imputed <- matrix(y, ns, nt)
      
      for (i in 1:ns) {
        n.thresh.miss <- sum(thresh.days[i, ]) + sum(missing.days[i, ])
        if (n.thresh.miss > 0){
          s.22     <- sig[-i, -i]
          s.22.inv <- chol2inv(chol(s.22))
          s.12     <- matrix(sig[i, -i], 1, (ns - 1))
          s.11     <- sig[i, i]
          for (t in 1:nt) {
          	e.y <- (mu[i, t]) + s.12 %*% s.22.inv %*% (y[-i, t] - mu[-i, t])
            s.y <- sqrt((s.11 - s.12 %*% s.22.inv $*% t(s.12)) / sigma[t])
            
            y.impute.t <- rTNorm(mn=e.y, sd=s.y, lower=-Inf, upper=thresh.mtx.k[i, ])
            
            # if any z come back as -Inf it's because P(Y < T) = 0                    
            if (y.impute.t == -Inf) {
              y.impute.t <- thresh.mtx.fudge
            }
              
            y.missing.t <- rnorm(n=1, mean=e.y, sd=s.y)
            
            if (thresh.days[i, t]) {
              y.imputed[i, t] <- y.impute.t
            }
            
            if (missing.days[i, t]) {
              y.imputed[i, t] <- y.missing.t
            }
          }  # end nt
        }  # fi n.thresh.miss > 0
      }  # end ns
      
      y <- y.imputed
    }  # fi thresh != 0
    
    # update random effects
    if (debug) {print("random z")}
    if (!fixz) {  # debug
      for (t in 1:nt) {
        for (k in 1:nknots) {
          these   <- which(partition[, t] == k)
          r1      <- y[these, t] - x.beta[these, t]
          r2      <- y[-these, t] - mu[-these, t]
          prec.11 <- prec[these, these]
          prec.21 <- prec[-these, these]
          mu.z    <- delta * sum(t(r1) * prec.11 + t(r2) * prec.21) / (1 - delta^2)
          prec.z  <- delta^2 * sum(prec.11) / (sigma * (1 - delta^2)) + 1 / sigma
          var.z   <- 1 / prec.z
        
          z.new <- rTNorm(mn=(var.z * mu.z), sd=sqrt(var.z), lower=0, upper=Inf) 
          z.knots[k, t] <- z.new
          z.sites[these, t] <- z.new
        }  # end nknots
      }  # end nt
    }  # fi !fixz
    
    # update partitions
    if (debug) { print("knots") }
    if (!fixknots) {  # debug
      if (nknots > 1) {
      	can.knots <- knots
        for (t in 1:nt) {
          att.w[t] <- att.w[t] + 1
          ss.t <- t(y[, t] - mu[, t]) %*% prec %*% (y[, t] - mu[, t])
          cur.ll <- -0.5 * ss.t / (sigma[t] * (1 - delta^2))
          
          for (k in 1:nknots) {
        	can.knots[[1]][k, ] <- rnorm(2, knots[[t]][k, ], mh.w[t])
          }
        
          can.partition <- Membership(s=s, knots=can.knots)
          can.z.sites <- ZBySites(z=z.knots, partition=can.partition)
          can.mu.t <- x.beta[, t] + delta * can.z.sites
          
          can.ss <- t(y[, t] - can.mu.t) %*% prec %*% (y[, t] - can.mu.t)
          can.ll <- -0.5 * can.ss / (sigma[t] * (1 - delta^2))
          
          rej <- sum(can.ll - cur.ll) # prior is uniform and candidate is symmetric
          
          if (!is.na(rej)) {if (-rej < rexp(1, 1)) {
            knots[[t]]     <- can.knots[[1]]
            partition[, t] <- can.partition
            z.sites[, t]   <- can.z.sites
            mu[, t]        <- can.mu.t
            acc.w[t]        <- acc.w[t] + 1
          }}
          
        } # end nt
        
      }  # fi nknots > 1
    }  # fi !fixknots
    
    # update beta
    if (debug) { print("beta") }
    if (!fixbeta) {  # debug
      prec.beta <- diag(p) / beta.s^2
      e.beta <- rep(beta.m, p)
      
      beta.post <- BetaPosterior(prec.beta=prec.beta, e.beta=e.beta, x=x, y=y, z.sites=z.sites,
                                 prec=prec, delta=delta, sigma=sigma, nt=nt)
                                 
      vvv <- beta.post$vvv
      mmm <- beta.post$mmm
      
      beta <- as.vector(vvv %*% mmm + t(chol(vvv)) %*% rnorm(n=p))
      
      for (t in 1:nt) {
  	    x.beta[, t] <- x[, t, ] %*% beta
  	    mu[, t] <- x.beta[, t] + delta * z.sites[, t]
      }
      
      
    }  #fi !fixbeta
    
    # update sigma
    if (debug) { print("sigma") }
    if (!fixsigma) {
      sigma.inv <- rep(0, nt)
      alpha.star <- sigma.a + nknots / 2 + ns / 2
      for (t in 1:nt) {
        beta.star <- sigma.b + sum(z[, t]^2) / 2 + 
                     t(y[, t] - mu[, t]) %*% prec %*% (y[, t] - mu[, t]) /  (2 * (1 - delta^2))
        sigma.inv[t] <- rgamma(n=1, shape=alpha.star, scale=beta.star)
      }
      
      sigma <- 1 / sigma.inv
     
    }  # fi !fixsigma
    
    # update delta
    if (debug) { print("delta") }
    if (!fixdelta) {
      cur.ll <- LLike(y=y, x.beta=x.beta, sigma=sigma, delta=delta, log.det=log.det, 
                     z.sites=z.sites, log=T)
      att.delta <- att.delta + 1
      
      alpha.skew <- delta / sqrt(1 - delta^2)
      can.alpha.skew <- rnorm(1, alpha.skew, mh.delta)
      can.delta <- can.alpha.skew / sqrt(1 + can.alpha.skew^2)
      
      can.ll <- LLike(y=y, x.beta=x.beta, sigma=sigma, delta=can.delta, 
                      log.det=log.det, z.sites=z.sites, log=T)
                     
      rej <- sum(can.ll - cur.ll)
      if (!is.na(rej)) { if(-rej < rexp(1, 1)) {
        delta <- can.delts
      }}
      
    }  # fi !fixdelta
    
  }  # end nthin
  }  # end iter   
  return(results)
}
