#########################################################################
# MCMC 
#
# TODO: Add in model description here
#
#ASSUMES x and y are in [0,1]^2
#########################################################################

mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL, 
                 thresh=0, thresh.quant=T, nknots=1, 
                 iters=5000, burn=1000, update=100, thin=1,
                 iterplot=F, plotname=NULL, method="t",
                 # initial values
                 beta.init=NULL, 
                 tau.init=1,
                 tau.alpha.init=0.1, 
                 tau.beta.init=0.1,
                 rho.init=0.5, 
                 nu.init=0.5, 
                 alpha.init=0.5,
                 # priors
                 beta.m=0, beta.s=10, 
                 tau.alpha.m=0, tau.alpha.s=1, 
                 tau.beta.a=0.1, tau.beta.b=0.1,
                 logrho.m=-2, logrho.s=1,
                 lognu.m=-1.2, lognu.s=1,
                 alpha.m=0, alpha.s=1,
                 # skew inits
                 z.init=1, z.alpha.init=0,
                 # skew priors
                 z.alpha.m=0, z.alpha.s=2, skew=T
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
  if (thresh.quant & thresh > 0) {  # threshold based on sample quantiles
    thresh.data <- quantile(y, thresh, na.rm=T)
  } else {
    thresh.data <- thresh
  }
  thresh.mtx  <- matrix(thresh.data, ns, nt)  # want as a matrix for easy replacement
  thresh.obs  <- !is.na(y) & (y < thresh.mtx)
  
  missing.obs <- is.na(y)
  y[missing.obs]  <- mean(y, na.rm=T)

    
  # initialize partition
  knots_con <- array(rnorm(nknots * nt * 2), c(nknots, 2, nt))
  knots     <- pnorm(knots_con)
    
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
    tau.init <- matrix(tau.init, ns, nt) 
  } else if (method == "t") {  # knots vary by day and partition
  	if (length(tau.init) == 1) {
  	  cat("\t initializing all tau terms to", tau.init, "\n")
  	}
  	tau  <- matrix(tau.init, nknots, nt)
    taug[, t] <- tau[g[, t], t]   
  }
  
  zg <- matrix(0, ns, nt)
  if (length(z.init) == 1) {
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
 
  C <- CorFx(d=d, alpha=alpha, rho=rho, nu=nu)
  prec.cor<-solve(C)
  
  # MH tuning params
  acc.w     <- att.w     <- mh.w     <- rep(0.1, nt)  # knot locations
  acc.tau   <- att.tau   <- mh.tau   <- matrix(1, nrow=nknots, ncol=nt)
  acc.delta <- att.delta <- mh.delta <- 0.1  
  acc.rho   <- att.rho   <- mh.rho   <- 0.1
  acc.nu    <- att.nu    <- mh.nu    <- 0.1
  acc.alpha <- att.alpha <- mh.alpha <- 0.5
  att.tau   <- acc.tau   <- mh.tau   <- rep(1,ns)    
  
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
  
  for (iter in 1:iters) { for (ttt in 1:thin) {
    
    # impute data below threshold
    if (thresh != 0) {
      mu <- x.beta + z.alpha * zg
      thresh.mtx.fudge <- 0.99999 * thresh.mtx  # numerical stability
      y.imputed <- matrix(y, ns, nt)
      
      cor <- alpha * matern(u=d, phi=rho, kappa=nu)
      for (t in 1:nt) {
        these.thresh.obs <- thresh.obs[, t]
        these.missing.obs <- missing.obs[, t]

        # spatial error
        sig.t <- 1 / sqrt(taug[, t])
        theta.t <- sig.t * t(chol(cor)) %*% rnorm(ns, 0, 1)  # generate for all sites
       
        # nugget error
        # new expected value and standard deviation
        e.y <- mu[, t] + theta.t
        s.y <- sqrt((1 - alpha) * sig.t)
        upper.y <- thresh.mtx[, t]
        
        y.impute.t <- rTNorm(mn=e.y, sd=s.y, lower=-Inf, upper=upper.y)
       
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
       g12    <- sqrt(taug[,t])
       prec.t <- diag(g12) %*% prec.cor %*% diag(g12) #faster with sweeps!
       ttt    <- t(x.t) %*% prec.t
       vvv    <- vvv + ttt %*% x.t
       mmm    <- mmm + ttt %*% (y[, t] - z.alpha * zg[, t])
    } 
    vvv <- solve(vvv)
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
        att.w[1]      <- att.w[1] + 1
        can_knots_con <- knots_con[, , t] + mh.w[1] * rnorm(2 * nknots)
        can_knots     <- pnorm(can_knots_con)
        cang          <- mem(s, can_knots)
        cantaug       <- tau[cang, t]
        canzg         <- z[cang, t]
        canres        <- y - x.beta - z.alpha * canzg   

        R <- -0.5 * quad.form(prec.cor, sqrt(cantaug) * canres[, t]) +
              0.5 * quad.form(prec.cor, sqrt(taug[, t]) * res[, t]) +
              0.5 * sum(log(cantaug)) -
              0.5 * sum(log(taug[, t])) +
              sum(dnorm(can_knots_con, log=TRUE)) -
              sum(dnorm(knots_con[, , t], log=TRUE))
          
        if (!is.na(R)) { if (log(runif(1)) < R) {
          knots_con[, , t] <- can_knots_con
          knots[, , t]     <- can_knots
          g[, t]           <- cang
          taug[, t]        <- cantaug
          zg[, t]          <- canzg
          acc.w[1]         <- acc.w[1] + 1
        }}
        
        # find the avg num of partitions per day
        for (k in 1:nknots) {
          if (sum(g[, t] == k) != 0) {
            avgparts[t] <- avgparts[t] + 1
          } 
        }
      }
    }
    
    

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

      if (nknots == 1) {
        for (t in 1:nt) {
          res.t <- res[, t]
          rss.t <- quad.form(prec.cor, res.t)
          tau[1, t] <- rgamma(1, tau.alpha + 0.5 * ns, tau.beta + 0.5 * rss.t)
          taug[, t] <- tau[1, t]
        }

      } else {

        for(t in 1:nt){
          res.t <- res[, t]
          curll[t] <- 0.5 * sum(log(taug[, t])) -
                      0.5 * quad.form(prec.cor, sqrt(taug[,t]) * res.t)

          for (k in 1:nknots) {

            these  <- which(g[, t] == k)
            nparts <- length(these)

            if(nparts==0){
              tau[k, t]<-rgamma(1, tau.alpha, tau.beta)
            }       
            else{      
              att.tau[nparts] <- att.tau[nparts] + 1
              aaa <- nparts / 2 + tau.alpha
              bbb <- quad.form(prec.cor[these, these], res.t[these]) / 2 + tau.beta

              aaa <- aaa / mh.tau[nparts]
              bbb <- bbb / mh.tau[nparts]

              cantau    <- tau[, t]
              cantau[k] <- rgamma(1, aaa, bbb)
              cantaug   <- cantau[g[, t]]

              canll <- 0.5 * sum(log(cantaug)) -
                       0.5 * quad.form(prec.cor, sqrt(cantaug) * res.t)

              R    <- canll - curll[t] +
                      dgamma(cantau[k], tau.alpha, tau.beta, log=TRUE) -
                      dgamma(tau[k, t], tau.alpha, tau.beta, log=TRUE) +
                      dgamma(tau[k, t], aaa, bbb, log=TRUE) -
                      dgamma(cantau[k], aaa, bbb, log=TRUE)

              if (!is.na(R)) { if (log(runif(1)) < R) {
                acc.tau[nparts] <- acc.tau[nparts] + 1
                tau[, t]  <- cantau
                taug[, t] <- cantaug
                curll[t]  <- canll
              }}
            }#end ifelse
          }#end k
        }  # end t
      }#end t
      
      # update tau.alpha and tau.beta
      a.star <- tau.beta.a + nt * nknots * tau.alpha
      b.star <- tau.beta.b + sum(tau)
      tau.beta <- rgamma(1, a.star, b.star)

      lll <- mmm <- seq(0.5, 10, 0.1)
      for (l in 1:length(lll)) {
        lll[l] <- sum(dgamma(tau, mmm[l], tau.beta, log=T))
      }
      tau.alpha <- sample(mmm, 1, prob=exp(lll - max(lll)))
    }
   
    
    # update rho and nu
    att.rho <- att.rho + 1
    att.nu  <- att.nu + 1
    
    logrho <- log(rho)
    can.logrho <- rnorm(1, logrho, mh.rho)
    can.rho <- exp(can.logrho)
    
    lognu  <- log(nu)
    can.lognu <- rnorm(1, lognu, mh.nu)
    can.nu <- exp(can.lognu)
    
    can.C <- CorFx(d=d, alpha=alpha, rho=can.rho, nu=can.nu)
    can.prec.cor <- solve(can.C)
    can.logdet.prec <- log(det(can.prec.cor))
    logdet.prec <- log(det(prec.cor))  
    
    can.rss <- rep(NA, nt)
    cur.rss <- rep(NA, nt)
    for (t in 1:nt) {
      can.rss[t] <- quad.form(can.prec.cor, sqrt(taug[, t]) * res[, t]) 
      cur.rss[t] <- quad.form(prec.cor, sqrt(taug[, t]) * res[, t])
    }
    
    R <- -0.5 * sum(can.rss - cur.rss) + 
          0.5 * nt * (can.logdet.prec - logdet.prec) + 
          dnorm(can.logrho, logrho.m, logrho.s, log=T) - 
          dnorm(logrho, logrho.m, logrho.s, log=T) + 
          dnorm(can.lognu, lognu.m, lognu.s, log=T) - 
          dnorm(lognu, lognu.m, lognu.s, log=T)
         
    if (!is.na(R)) { if (log(runif(1)) < R) {
      rho <- can.rho
      nu <- can.nu
      prec.cor <- can.prec.cor
      logdet.prec <- can.logdet.prec
      cur.rss <- can.rss
      acc.rho <- acc.rho + 1
      acc.nu  <- acc.nu + 1
    }}
    
    if (att.rho > 50) {
      if (acc.rho / att.rho < 0.25) { mh.rho <- mh.rho * 0.8 }
      if (acc.rho / att.rho > 0.50) { mh.rho <- mh.rho * 1.2 }
      acc.rho <- att.rho <- 1
    }
    
    if (att.nu > 50) {
      if (acc.nu / att.nu < 0.25) { mh.nu <- mh.nu * 0.8 }
      if (acc.nu / att.nu > 0.50) { mh.nu <- mh.nu * 1.2 }
      acc.nu <- att.nu <- 1
    }
       
    # update alpha
    att.alpha <- att.alpha + 1
    
    norm.alpha <- qnorm(alpha)
    can.norm.alpha <- rnorm(1, norm.alpha, mh.alpha)
    can.alpha <- pnorm(can.norm.alpha)
    
    can.C <- CorFx(d=d, alpha=can.alpha, rho=rho, nu=nu)
    can.prec.cor <- solve(can.C)
    can.logdet.prec <- log(det(can.prec.cor))
    logdet.prec <- log(det(prec.cor))
    
    can.rss <- rep(NA, nt)
    for (t in 1:nt) {
      can.rss[t] <- quad.form(can.prec.cor, sqrt(taug[, t]) * res[, t])
    }
    
    R <- -0.5 * sum(can.rss - cur.rss) + 
          0.5 * nt * (can.logdet.prec - logdet.prec) +
          dnorm(can.norm.alpha, mean=alpha.m, sd=alpha.s, log=T) - 
          dnorm(norm.alpha, mean=alpha.m, sd=alpha.s, log=T)
    
    if (!is.na(R)) { if (log(runif(1)) < R) {
      alpha <- can.alpha
      prec.cor <- can.prec.cor
      cur.rss <- can.rss
      acc.alpha <- acc.alpha + 1
    }}
    
    if (att.alpha > 50) {
      if (acc.alpha / att.alpha < 0.25) { mh.alpha <- mh.alpha * 0.8 }
      if (acc.alpha / att.alpha > 0.50) { mh.alpha <- mh.alpha * 1.2 }
      acc.alpha <- att.alpha <- 1
    }
    
    if (skew) {
      # Skewness parameter
      vvv <- 1 / z.alpha.s^2
      mmm <- z.alpha.m
    
      for (t in 1:nt) {
      	prec.cov <- diag(sqrt(taug[, t])) %*% prec.cor %*% diag(sqrt(taug[, t]))
        ttt <- zg[, t] %*% prec.cov
        vvv <- vvv + ttt %*% zg[, t]
        mmm <- mmm + ttt %*% (y[, t] - x.beta[, t])
      }
    
      vvv <- 1 / vvv
      mmm <- vvv * mmm 
    
      z.alpha <- rnorm(1, mmm, sqrt(vvv))
    
      # z random effect
      mu <- y - x.beta - z.alpha * zg
      
      for (t in 1:nt) {
        for (k in 1:nknots) {
          these <- which(g[, t] == k)
          r.1 <- y[these, t, drop=F] - x.beta[these, t, drop=F]
          r.2 <- y[-these, t, drop=F] - mu[-these, t, drop=F]
          prec.cov <- diag(sqrt(taug[, t])) %*% prec.cor %*% diag(sqrt(taug[, t]))
          prec.11 <- prec.cov[these, these, drop=F]  # with cov
          prec.21 <- prec.cov[-these, these, drop=F]
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
    }    
  }  # end nthin
  
  # predictions
  if (predictions) {
    s.11 <- as.matrix(alpha * matern(d11, rho, nu), np, np)
    diag(s.11) <- 1
    s.12 <- as.matrix(alpha * matern(d12, rho, nu), np, ns)
    s.12.22.inv <- s.12 %*% prec.cor
    corp <- s.11 - s.12.22.inv %*% t(s.12)
    corp <- diag(corp)
  
    yp <- matrix(NA, np, nt)

    for (t in 1:nt) {
      xp.beta  <- x.pred[, t, ] %*% beta
      if (nknots == 1) {
        gp <- 1
      } else {
        gp       <- mem(s.pred, knots[, , t])  # find the right partition
      }
      zgp    <- z[gp, t]
      siggp  <- 1 / sqrt(tau[gp, t])  # get the partition's variance
      taug.t <- sqrt(taug[, t])

      mup <- xp.beta + z.alpha * zgp - sqrt(siggp) * s.12.22.inv %*% (taug.t * res[, t])
      sdp <- sqrt(siggp * corp)
      
      yp[, t] <- mup + sdp * rnorm(np, 0, 1)
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
  
  if (iter %% update == 0) {
    if (iterplot) {
      if (skew) {
        par(mfrow=c(3, 3))
      } else {
        par(mfrow=c(3, 2))
      }
      plot(keepers.beta[1:iter, 1], type="l")
      plot(keepers.tau.alpha[1:iter], type="l")
      plot(keepers.tau.beta[1:iter], type="l")
      plot(keepers.rho[1:iter], type="l")
      plot(keepers.nu[1:iter], type="l")
      plot(keepers.alpha[1:iter], type="l")
      if (skew) {
        plot(keepers.z.alpha[1:iter], type="l")
        plot(keepers.z[1:iter, 1, 1], type="l")
        plot(keepers.z[1:iter, 1, 48], type="l")
      }
    }
    cat("\t iter", iter, "\n")
  }
  
  
  } #end iters

if (nknots == 1) {
  keepers.avgparts = matrix(1, iters, nt)
}
if (!predictions) {
  y.pred <- NULL
}
if (!skew) {
  z <- NULL
  z.alpha <- NULL
}

results <- list(tau=keepers.tau, 
                beta=keepers.beta,
                tau.alpha=keepers.tau.alpha,
                tau.beta=keepers.tau.beta,
                rho=keepers.rho,
                nu=keepers.nu,
                alpha=keepers.alpha,
                yp=y.pred,
                z.alpha=z.alpha,
                z=z,
                avgparts=keepers.avgparts)

return(results)
}#end mcmc()

