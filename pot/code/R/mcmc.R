#########################################################################
# MCMC 
#
# TODO: Add in model description here
#########################################################################

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL, 
                 thresh=0, thresh.quant=T, nknots=1,          
                 iters=5000, burn=1000, update=100, thin=1, scale=T,
                 iterplot=F, plotname=NULL,
                 # initial values
                 beta.init=NULL, sigma.init=1,
                 rho.init=0.5, nu.init=0.5, alpha.init=0.5,
                 delta.init=0,
                 # priors
                 beta.m=0, beta.s=10, sigma.a=0.1, sigma.b=10,
                 logrho.m=-2, logrho.s=1,
                 lognu.m=-1.2, lognu.s=1,
                 # debugging settings
                 debug=F, knots.init, z.init, 
                 fixknots=F, fixz=F, fixbeta=F, fixsigma=F,
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
  	z.sites[, t] <- ZBySites(z.knots[, t], partition[, t])
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
  
  if (length(sigma.init) == 1) {
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
  acc.delta <- att.delta <- mh.delta <- 10  
  acc.rho   <- att.rho   <- mh.rho   <- 0.1
  acc.nu    <- att.nu    <- mh.nu    <- 0.1
  acc.alpha <- att.alpha <- mh.alpha <- 0.1
  
  # storage
  if (debug) { print("storage") }
  keepers.z <- array(NA, dim=c(iters, nknots, nt))
  keepers.beta <- matrix(NA, nrow=iters, ncol=p)
  keepers.sigma <- keepers.ll <- matrix(NA, nrow=iters, ncol=nt)
  keepers.delta <- keepers.rho <- keepers.nu <- keepers.alpha <- rep(NA, iters)
  
  # initial values
  mu <- x.beta + delta * z.sites
  cur.ll <- LLike(y=y, x.beta=x.beta, sigma=sigma, delta=delta, prec=prec, 
                  log.det=log.det, z.sites=z.sites, log=T)
  
  for (iter in 1:iters) { for (ttt in 1:thin) {
    
    # impute data below threshold
    if (debug) { print("impute") }
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
            s.y <- sqrt((s.11 - s.12 %*% s.22.inv %*% t(s.12)) / sigma[t])
            
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
          mu.z    <- delta * sum(t(r1) %*% prec.11 + t(r2) %*% prec.21) / (1 - delta^2)
          prec.z  <- delta^2 * sum(prec.11) / (sigma[t] * (1 - delta^2)) + 1 / sigma[t]
          var.z   <- 1 / prec.z
        
          z.new <- rTNorm(mn=(var.z * mu.z), sd=sqrt(var.z), lower=0, upper=Inf) 
          z.knots[k, t] <- z.new
          z.sites[these, t] <- z.new
        }  # end nknots
      }  # end nt
      # update means
      mu <- x.beta + delta * z.sites
    }  # fi !fixz
    
    # update partitions
    if (debug) { print("knots") }
    if (!fixknots) {  # debug
      if (nknots > 1) {
        for (t in 1:nt) {
          can.knots <- knots
          att.w[t] <- att.w[t] + 1
          ss.t <- t(y[, t] - mu[, t]) %*% prec %*% (y[, t] - mu[, t])
          cur.ll <- -0.5 * ss.t / (sigma[t] * (1 - delta^2))
          
          for (k in 1:nknots) {
        	can.knots[[t]][k, ] <- rnorm(2, knots[[t]][k, ], mh.w[t])
          }
        
          can.partition <- Membership(s=s, knots=can.knots)
          can.z.sites <- ZBySites(z=z.knots, partition=can.partition)
          can.mu.t <- x.beta[, t] + delta * can.z.sites
          
          can.ss <- t(y[, t] - can.mu.t) %*% prec %*% (y[, t] - can.mu.t)
          can.ll <- -0.5 * can.ss / (sigma[t] * (1 - delta^2))
          
          rej <- sum(can.ll - cur.ll) # prior is uniform and candidate is symmetric
          
          if (!is.na(rej)) {if (-rej < rexp(1, 1)) {
            knots[[t]]     <- can.knots[[t]]
            partition[, t] <- can.partition
            z.sites[, t]   <- can.z.sites
            mu[, t]        <- can.mu.t
            acc.w[t]       <- acc.w[t] + 1
            cur.ll         <- can.ll
          }}
          
        } # end nt
        
      }  # fi nknots > 1
    }  # fi !fixknots
    
    # update beta
    if (debug) { print("beta") }
    if (!fixbeta) {  # debug
      prec.beta <- diag(p) / beta.s^2
      e.beta <- rep(beta.m, p)
      
      beta.post <- BetaPosterior(prec.beta=prec.beta, e.beta=e.beta, x=x, y=y, z=z.sites,
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
        beta.star <- sigma.b + sum(z.knots[, t]^2) / 2 + 
                     t(y[, t] - mu[, t]) %*% prec %*% (y[, t] - mu[, t]) /  (2 * (1 - delta^2))
        sigma.inv[t] <- rgamma(n=1, shape=alpha.star, rate=beta.star)
      }
      
      sigma <- 1 / sigma.inv
     
    }  # fi !fixsigma
    
    # update delta
    if (debug) { print("delta") }
    if (!fixdelta) {
      att.delta <- att.delta + 1
      
      cur.ll <- LLike(y=y, x.beta=x.beta, sigma=sigma, delta=delta, prec=prec,
                      log.det=log.det, z.sites=z.sites, log=T)
      
      alpha.skew <- delta / sqrt(1 - delta^2)
      can.alpha.skew <- rnorm(1, alpha.skew, mh.delta)
      can.delta <- can.alpha.skew / sqrt(1 + can.alpha.skew^2)
      
      can.ll <- LLike(y=y, x.beta=x.beta, sigma=sigma, delta=can.delta, 
                      prec=prec, log.det=log.det, z.sites=z.sites, log=T)
                     
      rej <- sum(can.ll - cur.ll)
      if (!is.na(rej)) { if(-rej < rexp(1, 1)) {
        delta  <- can.delta
        mu <- x.beta + delta * z.sites
        cur.ll <- can.ll
        acc.delta <- acc.delta + 1
      }}
      
    }  # fi !fixdelta
    
    # rho and nu
    if (debug) { print("rho and nu") }
    if (!fixrho || !fixnu) {

      if (!fixrho) {
      	att.rho    <- att.rho + 1
        logrho     <- log(rho)
        can.logrho <- rnorm(1, logrho, mh.rho)
      } else {
        can.logrho <- log(rho)
      }
      can.rho <- exp(can.logrho)
      
      if (!fixnu) {
      	att.nu    <- att.nu + 1
        lognu     <- log(nu)
        can.lognu <- rnorm(1, lognu, mh.nu)
      } else {
        can.lognu <- log(nu)
      }
      can.nu <- exp(can.lognu)
      
      if (can.nu <= 10) {  # for numerical stability
        can.cor.mtx <- SpatCor(d=d, alpha=alpha, rho=can.rho, nu=can.nu)
        can.sig     <- can.cor.mtx$sig
        can.prec    <- can.cor.mtx$prec
        can.log.det <- can.cor.mtx$log.det      
      
      
        can.ll <- LLike(y=y, x.beta=x.beta, sigma=sigma, delta=delta,
                        prec=can.prec, log.det=can.log.det, z.sites=z.sites, log=T)
                      
        rej <- sum(can.ll - cur.ll) + 
               dnorm(can.logrho, logrho.m, logrho.s, log=T) - 
               dnorm(logrho, logrho.m, logrho.s, log=T) +
               dnorm(can.lognu, lognu.m, lognu.s, log=T) -
               dnorm(lognu, lognu.m, lognu.s, log=T)
      
        if (!is.na(rej)) { if (-rej < rexp(1, 1)) {
          rho     <- can.rho
          nu      <- can.nu
          cor.mtx <- can.cor.mtx
          sig     <- can.sig
          prec    <- can.prec
          log.det <- can.log.det
          cur.ll  <- can.ll
          if (!fixrho) { acc.rho <- acc.rho + 1 }
          if (!fixnu) { acc.nu <- acc.nu + 1}
        }} 
      } # fi can.nu <= 10
    } # fi !fixrho || !fixnu
    
    # alpha
    if (debug) { print("alpha") }
    if (!fixalpha) {
      att.alpha      <- att.alpha + 1
      norm.alpha     <- qnorm(alpha)
      can.norm.alpha <- rnorm(1, norm.alpha, mh.alpha)
      can.alpha      <- pnorm(can.norm.alpha)
      
      can.cor.mtx <- SpatCor(d=d, alpha=can.alpha, rho=rho, nu=nu)
      can.sig     <- can.cor.mtx$sig
      can.prec    <- can.cor.mtx$prec
      can.log.det <- can.cor.mtx$log.det
      
      can.ll <- LLike(y=y, x.beta=x.beta, sigma=sigma, delta=delta,
                      prec=can.prec, log.det=can.log.det, z.sites=z.sites, log=T)
                      
      rej <- sum(can.ll - cur.ll) + 
             dnorm(can.norm.alpha, log=T) - dnorm(norm.alpha, log=T)
      
      if (!is.na(rej)) { if (-rej < rexp(1, 1)) {
        alpha     <- can.alpha
        cor.mtx   <- can.cor.mtx
        sig       <- can.sig
        prec      <- can.prec
        log.det   <- can.log.det
        cur.ll    <- can.ll
        acc.alpha <- acc.alpha + 1
      }}
      
    }  # fi !fixalpha
    
  }  # end nthin
  
  ##############################################
  # MH adapt candidate distributions
  ##############################################
  if (debug) { print("MH update") }
  if (iter < burn / 2) { 
    for (t in 1:nt) {
      if (att.w[t] > 50) {
        if (acc.w[t] / att.w[t] < 0.25) { mh.w[t] <- mh.w[t] * 0.8 }
        if (acc.w[t] / att.w[t] > 0.50) { mh.w[t] <- mh.w[t] * 1.2 }
      }
    }
    
    if (att.delta > 50) {
      if (acc.delta / att.delta < 0.25) { mh.delta <- mh.delta * 0.8 }
      if (acc.delta / att.delta > 0.50) { mh.delta <- mh.delta * 1.2 }
      if (acc.rho / att.rho < 0.25) { mh.rho <- mh.rho * 0.8 }
      if (acc.rho / att.rho > 0.50) { mh.rho <- mh.rho * 1.2 }
      if (acc.nu / att.nu < 0.25) { mh.nu <- mh.nu * 0.8 }
      if (acc.nu / att.nu > 0.50) { mh.nu <- mh.nu * 1.2 }
      if (acc.alpha / att.alpha < 0.25) { mh.alpha <- mh.alpha * 0.8 }
      if (acc.alpha / att.alpha > 0.50) { mh.alpha <- mh.alpha * 1.2 }
    }
  }  # fi iter < burn / 2
  
  ##############################################
  # Spatial Predictions
  ##############################################
  if (np > 0) {
    partition.pred <- Membership(s=s.pred, knots=knots)
    
    for (t in 1:nt) {
      x.beta.pred  <- x.pred[, t, ] %*% beta
      z.sites.pred <- ZBySites(z=z.knots[, t], partition.pred)
      mu.pred      <- x.beta.pred + delta * z.sites.pred
      
      s.11     <- 1
      s.12     <- matrix(alpha * matern(d12, phi=rho, kappa=nu), nrow=np, ncol=ns)
      s.22.inv <- prec
      
      e.y.pred <- mu.pred - s.12 %*% s.22.inv %*% (y[, t] - mu[, t])
      v.y.pred <- sigma[t] * (1 - delta^2) * (s.11 - s.12 %*% s.22.inv %*% t(s.12))
      
      if (np > 1) {
        v.y.pred <- diag(diag(v.y.pred))
      }
      
      y.pred[, t, iter] <- rmvnorm(1, mean=e.y.pred, sigma=v.y.pred)
      
    }
  }  # fi np > 0
  
  ##############################################
  # Keep track of iterations
  ##############################################
  if (debug) { print("keepers") }
  keepers.z[iter, , ]   <- z.knots
  keepers.beta[iter, ]  <- beta
  keepers.sigma[iter, ] <- sigma
  keepers.delta[iter]   <- delta
  keepers.rho[iter]     <- rho
  keepers.nu[iter]      <- nu
  keepers.alpha[iter]   <- alpha
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
  	  	
  	  par(mfrow=c(3, 4))
  	  if (iter > burn) {
  	    start <- burn
  	  } else {
  	    start <- 1
  	  }
  	  plot(keepers.beta[start:iter, 1], ylab="beta0", xlab="iteration", 
           type="l")
      plot(keepers.beta[start:iter, 2], ylab="beta1", xlab="iteration",
           type="l", main=plotmain)
      plot(keepers.beta[start:iter, 3], ylab="beta2", xlab="iteration",
           type="l")
      plot(keepers.ll[start:iter], ylab="loglike", xlab="iteration",
           type="l")
      plot(keepers.delta[start:iter], ylab="delta", xlab="iteration", 
           type="l", main=bquote("ACCR" == .(accrate.delta)))
      plot(keepers.rho[start:iter], ylab="rho", xlab="iteration", 
            type="l", main=bquote("ACCR" == .(accrate.rho)))
      plot(keepers.nu[start:iter], ylab="nu", xlab="iteration", 
           type="l", main=bquote("ACCR" == .(accrate.nu)))
      plot(keepers.alpha[start:iter], ylab="alpha", xlab="iteration", 
           type="l", main=bquote("ACCR" == .(accrate.alpha)))
      plot(keepers.sigma[start:iter, 1], ylab="sigma 1", xlab="iteration", 
           type="l")
      plot(keepers.sigma[start:iter, 3], ylab="sigma 3", xlab="itertion",
           type="l")
      plot(keepers.z[start:iter, 1, 1], ylab="z 11", xlab="iteration", 
           type="l")
      plot(keepers.z[start:iter, 1, 3], ylab="z 31", xlab="iteration",
           type="l")
  	}
  	
  	if ((iter == iters) && !is.null(plotname)) {
  	  dev.off()
  	}
  }  # fi iterplot
  
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
                    sigma=keepers.sigma,
                    delta=keepers.delta,
                    rho=keepers.rho,
                    nu=keepers.nu,
                    alpha=keepers.alpha,
                    yp=y.pred)
  
  return(results)
}
