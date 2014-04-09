#########################################################################
# MCMC 
#
# Theta = (rt, mu, sig, xi, rho, nu, alpha)
# Model: Yts ~MVN(mu.y, Sig)
#        Sig = Matern(sqrt(rt), rho, nu, alpha)
#        rt ~iid IG(xi.r, sig.r)
#
# Priors:
#        betas.y ~ N(0, 10)        
#        logrho.y ~ N(-2, 1)
#        lognu.y ~ N(-1, 1)
#        alpha.y ~ unif(0, 1)
#        logsig.r ~ N(0, 1)
#        logxi.r ~ N(0, 0.1)
#    
#########################################################################

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL, thresh=1,
                 r.model = "gamma", #also allow "fixed"
                 beta.y.m=0, beta.y.s=10, beta.init=0,
                 logrho.y.m=-2, logrho.y.s=1, logrho.y.init=log(0.1),
                 lognu.y.m=-1.2, lognu.y.s=1, lognu.y.init=log(0.5),
                 alpha.y.a=1, alpha.y.b=1, alpha.y.init=0.5,    
                 logsig.r.m=0, logsig.r.s=1, logsig.r.init=0,
                 logxi.r.m=0, logxi.r.s=0.1, logxi.r.init=0,
                 nknots=1,          
                 iters=5000, burn=1000, update=100, thin=1, scale=T,
                 debug=F, iterplot=F, plotname=NULL,
                 fixbeta=F, fixr=F, fixrho=F, fixnu=F, fixalpha=F,
                 fixxir=F, fixsigr=F, fixknots=F,
                 r.inv.init, knots.init){
    
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
  thresh.quant <- quantile(y, thresh)
  thresh.mtx   <- matrix(thresh.quant, ns, nt)  # need thresh in matrix form
  thresh.obs   <- !is.na(y) & (y <= thresh.mtx)
    
  # if there are some missing days, store the indices before
  # the first imputation. then set initial value as the average ozone
  # for all sites on the day.
  missing.obs <- is.na(y)
  for (t in 1:nt) {
    missing.sites       <- missing.obs[, t]
    y[missing.sites, t] <- mean(y[, t], na.rm=T)
  }
    
  # setup the initial covariance matrix
  # create partition
  knots.min.1 <- range( s[,1] )[1]
  knots.max.1 <- range( s[,1] )[2]
  knots.min.2 <- range( s[,2] )[1]
  knots.max.2 <- range( s[,2] )[2]
  knots       <- matrix(NA, nknots, 2)
  knots[,1]   <- runif(n=nknots, min=knots.min.1, max=knots.max.1)
  knots[,2]   <- runif(n=nknots, min=knots.min.2, max=knots.max.2)
    
  if(fixknots){knots <- matrix(knots.init, nknots, 2)}  # debug
       
  member     <- Membership(s=s, knots=knots, y=y, x=x)
  partition  <- member$partition
  y.by.knots <- member$y.by.knots
  x.by.knots <- member$x.by.knots
            
  # parameters
  beta.y <- rep(0, p)
  beta.y[1] <- mean(y)  # set initial intercept
  if(fixbeta){beta.y <- beta.init}  #debug
   
  mu.y <- ExpectY(partition=partition, x.by.knots=x.by.knots, 
                  beta.y=beta.y, nknots=nknots)
   	
  logrho.y <- logrho.y.init
  rho.y    <- exp(logrho.y)
  lognu.y  <- lognu.y.init
  nu.y     <- exp(lognu.y)
  alpha.y  <- alpha.y.init
            
  logsig.r <- logsig.r.init
  sig.r    <- exp(logsig.r.init)
  logxi.r  <- logxi.r.init
  xi.r     <- exp(logxi.r.init)
    
  # random effect
  r.inv <- matrix(rep(1, nknots*nt), nknots, nt)  # one r for each knot and day
    
  if(fixr){r.inv <- matrix(r.inv.init, nknots, nt)}   # debug

  # correlation matrix blocked by partition
  cor.mtx <- SpatCor(d=d, alpha=alpha.y, rho=rho.y, nu=nu.y,
    			   	 partition=partition, nknots=nknots)
  sig <- cor.mtx$sig
  prec <- cor.mtx$prec
  log.det <- cor.mtx$log.det

  # MH tuning parameters
  acck <- attk <- mh.k <- 0.1          # knot locations
  accy <- atty <- mh.y <- rep(0.1, 3)  # model params
  accr <- attr <- mh.r <- 0.1          # hyperparameters
  mh.y[1] <- 0.05
  mh.y[3] <- 0.5

  # storage
  # print("storage")  
  keepers.r     <- array(0, dim=c(iters, nknots, nt))  # random effects
  keepers.knots <- array(0, dim=c(iters, 2, nknots))   # knot locations
  keepers.y     <- matrix(0, iters, 6)                 # model params + llike
  keepers.beta  <- matrix(0, iters, p)                 # lm for mu.y
    
  # initial values
  ss <- SumSquares(y.by.knots=y.by.knots, mu.y=mu.y, prec=prec, 
                   partition=partition, nt=nt, nknots=nknots)
  curll <- LLike(ss=ss, log.det=log.det, r.inv=r.inv, partition=partition, log=T)   
  plotmain <- plotname
  plotname.file <- paste("plots/", plotname, sep="")
  
  for (iter in 1:iters) { for(ttt in 1:thin) { 
    if (debug) { print("knots") }
    if (!fixknots) { #debug
    if (nknots > 1) {
      attk <- attk + 1
      canknots <- matrix(NA, nknots, 2)
      for (k in 1:nknots) {
        canknots[k, ] <- rnorm(2, knots[k, ], mh.k)
      }
      
      canmember     <- Membership(s=s, knots=canknots, y=y, x=x)
      canpartition  <- canmember$partition
      cany.by.knots <- canmember$y.by.knots
      canx.by.knots <- canmember$x.by.knots
      canmu.y       <- ExpectY(partition=canpartition, x.by.knots=canx.by.knots, 
                               beta.y=beta.y, nknots=nknots)
      cancor.mtx    <- SpatCor(d=d, alpha=alpha.y, rho=rho.y, nu=nu.y, 
                               partition=canpartition, nknots=nknots)
      canprec       <- cancor.mtx$prec
      canlog.det    <- cancor.mtx$log.det
      canss         <- SumSquares(y.by.knots=cany.by.knots, mu.y=canmu.y, 
                                  prec=canprec, partition=canpartition, 
                                  nt=nt, nknots=nknots)
      canll         <- LLike(ss=canss, log.det=canlog.det, r.inv=r.inv, 
                             partition=canpartition, log=T)
    		
      if(debug) { print(canknots) }
      
      rej <- sum(canll - curll) # prior is uniform and candidate is symmetric
      if (!is.na(rej)) {if (-rej < rexp(1, 1)) {
        knots      <- canknots
        member     <- canmember
        partition  <- canpartition
        y.by.knots <- cany.by.knots
        x.by.knots <- canx.by.knots
        mu.y       <- canmu.y     
        cor.mtx    <- cancor.mtx
        prec       <- canprec
        log.det    <- canlog.det
        sig        <- cancor.mtx$sig
        ss         <- canss
        curll      <- canll
        acck       <- acck + 1
      } }
		        
    } 
    } #fixknots     
        
		
    # impute missing and data below threshold
    if(debug) { print("impute") }     
    if (thresh != 0) { 
        	
      for (k in 1:nknots) {
        these      <- which(partition == k)
        ns.k       <- length(these)
        sig.k      <- sig[[k]]
        y.k        <- y.by.knots[[k]]
        mu.y.k     <- mu.y[[k]]
        r.inv.k    <- r.inv[k, ]
        thresh.mtx.k <- thresh.mtx[these, ]
         	
    	thresh.mtx.fudge <- 0.99999 * thresh.mtx.k
    	
        thresh.days.k  <- matrix(thresh.obs[these, ], ns.k, nt)
        missing.days.k <- matrix(missing.obs[these, ], ns.k, nt)
        y.k.imputed <- matrix(y.k, ns.k, nt)
        
        if (ns.k == 1){
          n.thresh.miss <- sum(thresh.days.k[1, ]) + sum(missing.days.k[1, ])
          if (n.thresh.miss > 0) {
            for (t in 1:nt) {
          	
          	  e.y <- mu.y.k[, t]
              s.y <- sqrt(sig.k / r.inv.k[t])
          
              y.impute.t  <- rTNorm(mn=e.y, sd=s.y, lower=-Inf,
                                    upper=thresh.mtx.k)
            
              # if any z come back as -Inf it's because P(Y < T) = 0
              if (y.impute.t == -Inf) {
                y.impute.t <- thresh.mtx.fudge
              }
                        
              y.missing.t <- rnorm(n=1, mean=e.y, sd=s.y)
          
              if (thresh.days.k[1, t]) {
                y.k.imputed[1, t] <- y.impute.t
              }
            
              if (missing.days.k[1, t]) {
                y.k.imputed[1, t] <- y.missing.t
              }
            }
          }  # end if n.thresh.miss > 0
        } else if (ns.k > 1) {
          for (i in 1:ns.k) {
            n.thresh.miss <- sum(thresh.days.k[i, ]) + sum(missing.days.k[i, ])
            if (n.thresh.miss > 0){
            s.22.k     <- sig.k[-i, -i]
            s.22.inv.k <- chol2inv(chol(s.22.k))
            s.12.k     <- matrix(sig.k[i, -i], 1, (ns.k - 1))
            s.11.k     <- sig.k[i, i]
              for (t in 1:nt) {
		        # s.22.inv.kt <- s.22.inv * r.inv.k[t]
		        # s.12.kt     <- s.12 / r.inv.k[t]
		        # s.11.kt     <- s.11 / r.inv.k[t]
		        # MVN impute 
                e.y      <- mu.y.k[i, t] + s.12.k %*% s.22.inv.k %*% 
                            (y.k[-i, t] - mu.y.k[-i, t])

                s.y      <- sqrt((s.11.k - s.12.k %*% s.22.inv.k %*% t(s.12.k))/r.inv.k[t])
		            
                y.impute.t  <- rTNorm(mn=e.y, sd=s.y, lower=-Inf,
                                      upper=thresh.mtx.k[i,])
              
                # if any z come back as -Inf it's because P(Y < T) = 0                    
                if (y.impute.t == -Inf) {
                  y.impute.t <- thresh.mtx.fudge
                }
              
                y.missing.t <- rnorm(n=1, mean=e.y, sd=s.y)
		                 
                # y.impute  <- Z2Y(z.impute, mu.y.k[i, ], r.inv.k)
                # y.missing <- Z2Y(z.missing, mu.y.k[i, ], r.inv.k)
              
                if (thresh.days.k[i, t]) {
                  y.k.imputed[i, t] <- y.impute.t
                }
              
                if (missing.days.k[i, t]) {
                  y.k.imputed[i, t] <- y.missing.t
                }
              
                #y.k[i, thresh.days]  <- y.impute[thresh.days]
                #y.k[i, missing.days] <- y.missing[missing.days]
		              
              }  # for t in 1:nt
            }  # end n.thresh.miss > 0
          }  # for i in 1:ns.k
        }  # if ns.k > 0
        y.by.knots[[k]] <- y.k.imputed
      } # end k in 1:nknots
    }  # end if thresh != 0
                
    ##############################################
    # Store updated values
    ##############################################
    mu.y  <- ExpectY(partition=partition, x.by.knots=x.by.knots, beta.y=beta.y,
                     nknots=nknots)
    ss    <- SumSquares(y.by.knots=y.by.knots, mu.y=mu.y, prec=prec,
                        partition=partition, nt=nt, nknots=nknots)
    curll <- LLike(ss=ss, log.det=log.det, r.inv=r.inv, partition=partition, log=T)
        
    ##########################################
    # random effects
    ##########################################
    if(debug) { print("gamre") }  
    
    if (r.model == "fixed") {
      r.inv <- rgamma(nknots, ns * nt / 2 + 0.1, sum(ss) / 2 + 0.1)
      r.inv <- matrix(r.inv, nrow=nknots, ncol=nt)
    } else if (r.model == "gamma") {
      # Sample r
      if (!fixr) { #debug
      r.inv <- matrix(NA, nrow=nknots, ncol=nt) 
      for (k in 1:nknots) {
        these      <- which(partition == k)
        ns.k       <- length(these)
        r.inv[k, ] <- rgamma(nt, ns.k / 2 + xi.r, ss[k, ] / 2 + sig.r)
      }
      }
            
      if(!fixsigr){ #debug
      # Update hyperparameters
      sig.r <- rgamma(1, xi.r * nt * nknots + 0.1, 0.1 + sum(r.inv))
      } #debug
			
      if(!fixxir){ #debug
      attr       <- attr + 1
      logxi.r    <- log(xi.r)
      canlogxi.r <- rnorm(1, logxi.r, mh.r)
      canxi.r    <- exp(canlogxi.r)
            
      curll.r <- sum(dgamma(r.inv, xi.r, sig.r, log=T))
      canll.r <- sum(dgamma(r.inv, canxi.r, sig.r, log=T))
			
      rej <- sum(canll.r - curll.r) +
             dnorm(canlogxi.r, logxi.r.m, logxi.r.s, log=T) -
             dnorm(logxi.r, logxi.r.m, logxi.r.s, log=T)
			
      if (!is.na(rej)) {if (-rej < rexp(1, 1)) {
        logxi.r <- canlogxi.r
        xi.r    <- canxi.r	
        accr    <- accr + 1
      } }
			
      } #debug
    } # end r.model==gamma
        
    curll <- LLike(ss=ss, log.det=log.det, r.inv=r.inv, partition=partition, log=T)        
        
    ##########################################
    # Expected value
    ##########################################
    # mu.y 
    if(debug) { print("mu.y") } 
    if (!fixbeta) { #debug
    
    prec.beta <- diag(p) / (beta.y.s^2)
    e.beta   <- rep(beta.y.m, p)
        
    beta.post <- BetaPosterior(prec.beta=prec.beta, e.beta=e.beta, 
                               x.by.knots=x.by.knots, y.by.knots=y.by.knots,
                               prec=prec, partition=partition, r.inv=r.inv, nt=nt)
    vvv <- beta.post$vvv
    mmm <- beta.post$mmm 

    beta.y <- vvv %*% mmm + t(chol(vvv)) %*% rnorm(p)
    beta.y <- as.vector(beta.y) 
    
    mu.y  <- ExpectY(partition=partition, x.by.knots=x.by.knots, beta.y=beta.y,
                     nknots=nknots)
    ss    <- SumSquares(y.by.knots=y.by.knots, mu.y=mu.y, prec=prec,
                        partition=partition, nt=nt, nknots=nknots)
    curll <- LLike(ss=ss, log.det=log.det, r.inv=r.inv, partition=partition, log=T)
       
    } # debug
		
	##########################################
    # Spatial covariance parameters
    ##########################################
        
    # rho and nu
    if(debug) { print("logrho.y and lognu.y") }
    atty[1] <- atty[1] + 1
    atty[2] <- atty[2] + 1
    
    if (!fixrho) {    
    	canlogrho.y <- rnorm(1, logrho.y, mh.y[1])
    } else {
    	canlogrho.y <- logrho.y
    }
    canrho.y    <- exp(canlogrho.y)
    
    if (!fixnu) {
    	canlognu.y <- rnorm(1, lognu.y, mh.y[2])
    } else {
    	canlognu.y <- lognu.y
    }
    cannu.y     <- exp(canlognu.y) 
    
    if(cannu.y <= 10){ # numerical stability
    cancor.mtx  <- SpatCor(d=d, alpha=alpha.y, rho=canrho.y, nu=cannu.y, 
                           partition=partition, nknots=nknots)
    canprec     <- cancor.mtx$prec
    canlog.det  <- cancor.mtx$log.det
    canss       <- SumSquares(y.by.knots=y.by.knots, mu.y = mu.y, 
                              prec=canprec, partition=partition, nt=nt, nknots=nknots)
    canll       <- LLike(ss=canss, log.det=canlog.det, r.inv=r.inv, 
                         partition=partition, log=T)
        
    rej <- sum(canll - curll) +
           dnorm(canlogrho.y, logrho.y.m, logrho.y.s, log=T) -
           dnorm(logrho.y, logrho.y.m, logrho.y.s, log=T) +
           dnorm(canlognu.y, lognu.y.m, lognu.y.s, log=T) -
           dnorm(lognu.y, lognu.y.m, lognu.y.s, log=T)
        
    if (!is.na(rej)) {if (-rej < rexp(1,1)) {
      logrho.y <- canlogrho.y
      rho.y    <- canrho.y
      lognu.y  <- canlognu.y
      nu.y     <- cannu.y
      cor.mtx  <- cancor.mtx
      prec     <- canprec
      log.det  <- canlog.det
      sig      <- cancor.mtx$sig
      ss       <- canss
      curll    <- canll
      accy[1]  <- accy[1] + 1
      accy[2]  <- accy[2] + 1
    }}
    }
        
    # alpha
    if(debug) { print("alpha.y") }
    if(!fixalpha){
      atty[3]    <- atty[3] + 1
      temp       <- qnorm(alpha.y)
      cantemp    <- rnorm(1, temp, mh.y[3])
      canalpha.y <- pnorm(cantemp)
      
      cancor.mtx <- SpatCor(d=d, alpha=canalpha.y, rho=rho.y, nu=nu.y, 
                            partition=partition, nknots=nknots)
      canprec    <- cancor.mtx$prec
      canlog.det <- cancor.mtx$log.det
      canss      <- SumSquares(y.by.knots=y.by.knots, mu.y=mu.y,
                               prec=canprec, partition=partition, nt=nt, nknots=nknots)
      canll      <- LLike(ss=canss, log.det=canlog.det, r.inv=r.inv, 
                          partition=partition, log=T)
        
      rej <- sum(canll - curll) +
             dnorm(cantemp, log=T) - dnorm(temp, log=T)
        
      if(!is.na(rej)){if(-rej < rexp(1,1)){
          alpha.y <- canalpha.y
          cor.mtx <- cancor.mtx
          prec    <- canprec
          log.det <- canlog.det
          sig     <- cancor.mtx$sig
          ss      <- canss
          curll   <- canll
          accy[3] <- accy[3] + 1
      }}
    }
        
    } # end nthin
    
    ##############################################
    # MH adapt candidate distributions
    ##############################################
    if(debug){ print("MH update") }
    if ((iter < burn / 2) & (attk >50)) {
      if (acck/attk < 0.25) {mh.k <- mh.k * 0.8}
      if (acck/attk > 0.5) {mh.k <- mh.k * 1.2}
      acck <- attk <- 0
    }
  
    if ((iter < burn / 2) & (attr > 50)) {
      if (accr / attr < 0.25) {mh.r <- mh.r * 0.8}
      if (accr / attr > 0.5) {mh.r <- mh.r * 1.2}
      accr <- attr <- 0
    }
  
    for (j in 1:length(atty)) {
      if(iter < burn/2 & atty[j] > 50){
        if(accy[j]/atty[j]<0.25){mh.y[j] <- mh.y[j]*0.8}
        if(accy[j]/atty[j]>0.5){mh.y[j] <- mh.y[j]*1.2}
        accy[j] <- atty[j] <- 0
      }
    }
    
    ##############################################
    # Keep track of iterations
    ##############################################
    if(debug){ print("keepers") }
    keepers.r[iter, , ]     <- 1/r.inv
    keepers.knots[iter, , ] <- knots
    keepers.beta[iter, ]    <- beta.y
    keepers.y[iter, ]       <- c(rho.y, nu.y, alpha.y, sig.r, xi.r, 
                                 sum(curll))
    
    ##############################################
    # Spatial Predictions
    ##############################################
    if (np > 0) {
      member.pred     <- Membership(s=s.pred, knots=knots, y=y.pred[, , iter], x=x.pred)
      partition.pred  <- member.pred$partition
      x.pred.by.knots <- member.pred$x.by.knots
      mu.pred         <- ExpectY(partition=partition.pred, 
                                 x.by.knots=x.pred.by.knots, beta.y=beta.y, 
                                 nknots=nknots)
      if(debug){ print("predict") }
      
      for (k in 1:nknots) {
      	these.pred <- which(partition.pred == k)
      	n.pred.k   <- length(these.pred)
      	if (n.pred.k > 0) { # don't need to worry about partitions with no preds
      	  mu.y.k     <- mu.y[[k]]
      	  mu.pred.k  <- mu.pred[[k]]
      	  these.obs  <- which(partition == k)
          n.obs.k    <- length(these.obs)
          y.obs.k    <- y.by.knots[[k]]
          d12.k      <- d12[these.pred, these.obs]
          
          d12.k    <- matrix(d12.k, n.pred.k, n.obs.k)
                   
          r.inv.k    <- r.inv[k, ]
          
          for(t in 1:nt){
	        r.inv.kt <- r.inv.k[t]
	        s.22.inv <- r.inv.kt * prec[[k]]
	        s.11     <- 1 / r.inv.kt	
	         
	        if (n.obs.k > 0) {
	          s.12 <- matern(d12.k, phi=rho.y, kappa=nu.y)
	          s.12 <- alpha.y * s.12 / r.inv.kt
	          s.12 <- matrix(s.12, nrow=n.pred.k, ncol=n.obs.k)  
	          
	          e.y.pred <- mu.pred.k[, t] - s.12 %*% s.22.inv %*% 
	                      (y.obs.k[,t] - mu.y.k[, t])
	          v.y.pred <- s.11 - s.12 %*% s.22.inv %*% t(s.12)
	          
	          if (n.pred.k > 1) {	          
	            v.y.pred <- diag(diag(v.y.pred))
	          }
	        
	        } else {
	          e.y.pred <- rep(0, n.pred.k)
	          v.y.pred <- diag(s.11, n.pred.k)          	
	        }
	       
	        y.pred.kt <- rmvnorm(1, mean=e.y.pred, sigma=v.y.pred)
           	  
      	    y.pred[these.pred, t, iter] <- y.pred.kt
	      }
      	}  # end n.pred.k>0
      }  # end k in 1:nknots
    }
    
    ##############################################
    # Display current value
    ##############################################
    if(debug){ print("plot") }
    
    if (iterplot) {
      # if ((iter %% update) == 0) {
      if (iter == iters) {
        pdf(file=plotname.file)
        product <- keepers.y[1:iter, 1] * keepers.y[1:iter, 3]
        accrate.y    <- round(accy / atty, 3)
        accrate.xi.r <- round(accr / attr, 3)
        par(mfrow=c(3, 3))
        plot(keepers.beta[1:iter, 1], ylab="beta0", xlab="iteration", 
             type="l")
        plot(keepers.beta[1:iter, 2], ylab="beta1", xlab="iteration",
             type="l", main=plotmain)
        plot(keepers.beta[1:iter, 3], ylab="beta3", xlab="iteration",
             type="l")
        plot(keepers.y[1:iter, 1], ylab="rho.y", xlab="iteration", 
             type="l", main=bquote("ACCR" == .(accrate.y[1])))
        plot(keepers.y[1:iter, 2], ylab="nu.y", xlab="iteration", 
              type="l", main=bquote("ACCR" == .(accrate.y[2])))
        plot(keepers.y[1:iter, 3], ylab="alpha.y", xlab="iteration", 
             type="l", main=bquote("ACCR" == .(accrate.y[3])))
        plot(keepers.y[1:iter, 4], ylab="sig.r", xlab="iteration", 
             type="l")
        plot(keepers.y[1:iter, 5], ylab="xi.r", xlab="iteration", 
             type="l", main=bquote("ACCR" == .(accrate.xi.r)))
        plot(keepers.y[1:iter, 6], ylab="loglike", xlab="itertion",
             type="l", main="Log Likelihood")
        dev.off()
      } 
    }
    if ((iter %% update) == 0) {
      cat("   Iter", iter, "complete. \n")
    }
    if (debug) { print(paste("iter: ", iter)) }
  
  } # end iters
    
  ##############################################
  # Return output
  ##############################################
  stop.time <- proc.time()
  results   <- list(time=stop.time-start.time,
                    gpdre=keepers.r,
                    params=keepers.y,
                    knots=keepers.knots,
                    beta=keepers.beta,
                    yp=y.pred)
     
  return(results)
}
