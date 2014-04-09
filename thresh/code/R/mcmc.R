# Model
#
# logit[P(y>thresh)] = Xprob%*%beta1
# y|y>thresh ~ GPD(lower bound = thresh, 
#                  scale = exp(Xscale%*%beta2), 
#                  shape = Xshape%*%beta3)
#
# Make sure all Xs are matrices, not vectors!

Bayes_GPD <- function(y, X.prob, X.sig, X.xi, 
                      Xp.prob=NULL, Xp.sig=NULL, Xp.xi=NULL, 
                      prob.m=0, prob.s=10, fix.prob=F,
                      sig.m=-2, sig.s=1, fix.sig=F,
                      xi.m=0, xi.s=0.2, fix.xi=F,
                      prob.init=NULL, sig.init = NULL, xi.init=NULL,
                      fix.thresh = T, thresh.a = 0.5, thresh.b = 0.95,
                      thresh, iters=2500, burn=100, update=100, thin=10,
                      iterplot=F, debug=F) {

  p1 <- ncol(X.prob) # beta 1 length
  p2 <- ncol(X.sig)  # beta 2 length
  p3 <- ncol(X.xi)   # beta 3 length

  z <- ifelse (y > thresh, 1, 0)  # over threshold
  Z <- X.prob                     # all covariates for P(Y > T)

  y <- y - thresh       # set thresh = 0
  n.exceed <- sum(y > 0)
  X1 <- matrix(X.prob[y > 0, ], nrow=n.exceed, ncol=p1) # covariates (prob)
  X2 <- matrix(X.sig[y > 0, ], nrow=n.exceed, ncol=p2)  # covariates (sig)
  X3 <- matrix(X.xi[y > 0, ], nrow=n.exceed, ncol=p3)   # covariates (xi)
  y <- y[y > 0]         # observations above threshold

  beta.prob <- rep(0, p1)  # for P(Y > T)
  beta.sig  <- rep(0, p2)  # GPD scale
  beta.xi   <- rep(0, p3)  # GPD shape
  
  # initial values
  if (!is.null(prob.init)) {
  	beta.prob[1] <- logit(prob.init)
  } else {
    beta.prob[1] <- logit(mean(z))
  }
  
  if (!is.null(sig.init)) {
    beta.sig[1] <- log(sig.init)
  } else {
    beta.sig[1] <- log(max(y) / 2)
  }
  
  if (!is.null(xi.init)) {
    beta.xi[1] <- xi.init
  } else {
    beta.xi[1]   <- 0.1
  }
  
  # storage
  gelf.ghosh <- rep(0, iters)
  CRPS <- rep(0, iters)
  keep.beta.prob <- matrix(0, iters, p1)
  keep.beta.sig  <- matrix(0, iters, p2)
  keep.beta.xi   <- matrix(0, iters, p3)
  colnames(keep.beta.prob) <- colnames(X.prob)
  colnames(keep.beta.sig)  <- colnames(X.sig)
  colnames(keep.beta.xi)   <- colnames(X.xi)
  
  # predictions
  np <- 0
  predictions <- !is.null(Xp.prob) & !is.null(Xp.sig) & !is.null(Xp.xi)
  y.pred <- NULL
  if (predictions) {
  	np <- nrow(Xp.prob)
  	y.pred <- matrix(0, nrow=iters, np)
  }

  MH.prob <- acc.prob <- att.prob <- rep(0.1, p1)
  MH.sig  <- acc.sig  <- att.sig  <- rep(0.1, p2)
  MH.xi   <- acc.xi   <- att.xi   <- rep(0.1, p3)

  for (i in 1:iters) { for (ttt in 1:thin) {
  if (debug) {cat("MH.sig is ", MH.sig, "\n")}
    # P(Y > T)
    cur.prob <- expit(Z %*% beta.prob)  # probability of exceeding
    curl.prob <- sum(dbinom(z, 1, cur.prob, log=T))  # only need above/below
    
    if (!fix.prob){
      for (j in 1:p1) {
      
        att.prob[j] <- att.prob[j] + 1
        can.beta.prob    <- beta.prob
        can.beta.prob[j] <- rnorm(1, beta.prob[j], MH.prob[j])
        can.prob <- expit(Z %*% can.beta.prob)
        canl.prob <- sum(dbinom(z, 1, can.prob, log=T))

        R <- canl.prob - curl.prob +
             dnorm(can.beta.prob[j], prob.m, prob.s, log=T) -
             dnorm(beta.prob[j], prob.m, prob.s, log=T)
          
        if (!is.na(R)) { if (runif(1) < exp(R)) {
          beta.prob <- can.beta.prob
          cur.prob <- can.prob
          curl.prob <- canl.prob
          acc.prob[j] <- acc.prob[j] + 1
        } }
      }
    }
    
    cur.prob.exceed <- expit(X1 %*% beta.prob)  # need for GP function
    
    # GPD parameters
    if(debug){cat("beta.sig is ", beta.sig, "\n")}
    cur.sig  <- exp(as.vector(X2 %*% beta.sig))
    cur.xi   <- as.vector(X3 %*% beta.xi)
    curl.gpd <- dGP(y=y, scale=cur.sig, shape=cur.xi)
    
    # scale
    if (!fix.sig) {
      att.sig <- att.sig + 1
      can.beta.sig <- rnorm(p2, beta.sig, MH.sig)
      can.sig <- exp(as.vector(X2 %*% can.beta.sig))
      canl.gpd <- dGP(y=y, scale=can.sig, shape=cur.xi)
    
      R <- canl.gpd - curl.gpd +
           sum(dnorm(can.beta.sig, sig.m, sig.s, log=T) -
               dnorm(beta.sig, sig.m, sig.s, log=T))
    
      if (!is.na(R)) { if (runif(1) < exp(R)) {
        beta.sig <- can.beta.sig
        cur.sig <- can.sig
        curl.gpd <- canl.gpd
        acc.sig <- acc.sig + 1
      }  }
    }
    
    # for (j in 1:p2) {   
      # att.sig[j] <- att.sig[j] + 1
      # can.beta.sig    <- beta.sig
      # can.beta.sig[j] <- rnorm(1, beta.sig[j], MH.sig[j])
      # if(debug){cat("scale beta=", can.beta.sig, "\n")}
      # can.sig  <- exp(as.vector(X2 %*% can.beta.sig))
      # if(debug){cat("scale can=", can.sig, "\n")}
      # canl.gpd <- GP(y, can.sig, cur.xi)

      # R <- canl.gpd - curl.gpd +
           # dnorm(can.beta.sig[j], sig.m, sig.s, log=T) -
           # dnorm(beta.sig[j], sig.m, sig.s, log=T)
           
      # if(debug){cat("scale R=", R, "\n")}
          
      # if (!is.na(R)) { if (runif(1) < exp(R)) {
        # beta.sig <- can.beta.sig
        # cur.sig  <- can.sig
        # curl.gpd <- canl.gpd
        # acc.sig[j] <- acc.sig[j] + 1
      # } }
    # }


    # xi
    if (!fix.xi) {
      att.xi <- att.xi + 1
      can.beta.xi <- rnorm(p3, beta.xi, MH.xi)
      can.xi <- as.vector(X3 %*% can.beta.xi)
      canl.gpd <- dGP(y=y, scale=cur.sig, shape=can.xi)
    
      R <- canl.gpd - curl.gpd +
           sum(dnorm(can.beta.xi, xi.m, xi.s, log=T) -
               dnorm(beta.xi, xi.m, xi.s, log=T))
    
      if (!is.na(R)) { if (runif(1) < exp(R)) {
        beta.xi <- can.beta.xi
        cur.xi <- can.xi
        curl.gpd <- canl.gpd
        acc.xi <- acc.xi + 1
      }}
    }
    # for(j in 1:p3){
      # att.xi[j] <- att.xi[j] + 1
      # can.beta.xi <- beta.xi
      # can.beta.xi[j] <- rnorm(1, beta.xi[j], MH.xi[j])
      # if(debug){cat("xi beta=", can.beta.xi, "\n")}
      # can.xi <- as.vector(X3 %*% can.beta.xi)
      # if(debug){cat("xi can=", can.xi, "\n")}
      # canl.gpd <- GP(y, cur.sig, can.xi)

      # R <- canl.gpd - curl.gpd +
           # dnorm(can.beta.xi[j], xi.m, xi.s, log=T) -
           # dnorm(beta.xi[j], xi.m, xi.s, log=T)
      
      # if(debug){cat("xi R=", R, "\n")}
           
      # if (!is.na(R)) { if (runif(1) < exp(R)) {
        # beta.xi <- can.beta.xi
        # cur.xi <- can.xi
        # curl.gpd <- canl.gpd
        # acc.xi[j] <- acc.xi[j] + 1
      # } }
      
    # }   
  
  }#end thin


  #Tune the algorithm:
  if (i < burn) { 
    for (j in 1:p1) { if (att.prob[j] > 50) {
      if (acc.prob[j] / att.prob[j] < 0.4) { MH.prob[j] <- MH.prob[j] * 0.8 }
      if (acc.prob[j] / att.prob[j] > 0.6) { MH.prob[j] <- MH.prob[j] * 1.2 }
      acc.prob[j] <- att.prob[j] <- 0
    } }         
     
    for (j in 1:p2) { if (att.sig[j] > 50) {
      if (acc.sig[j] / att.sig[j] < 0.4) {MH.sig[j] <- MH.sig[j] * 0.8}
      if (acc.sig[j] / att.sig[j] > 0.6) {MH.sig[j] <- MH.sig[j] * 1.2}
      acc.sig[j] <- att.sig[j] <- 0
    }  }      
     
    for (j in 1:p3) { if (att.xi[j] > 50) {
      if (acc.xi[j] / att.xi[j] < 0.4) {MH.xi[j] <- MH.xi[j] * 0.8}
      if (acc.xi[j] / att.xi[j] > 0.6) {MH.xi[j] <- MH.xi[j] * 1.2}
      acc.xi[j] <- att.xi[j] <- 0
    } }
  }        

  # model fit
  # gelf.ghosh
  # CRPS

  keep.beta.prob[i, ] <- beta.prob
  keep.beta.sig[i, ]  <- beta.sig
  keep.beta.xi[i, ]   <- beta.xi
  
  if (debug) {cat("mean y is ", mean(y) )}
  
  #################################################:
  ##########      Prediction          #############:
  #################################################:

  if (np > 0) {
  	cur.prob.p <- expit(Xp.prob %*% beta.prob)
    cur.sig.p <- exp(as.vector(Xp.sig %*% beta.sig))
    cur.xi.p <- as.vector(Xp.xi %*% beta.xi)
  	
  	u.p <- runif(np, 0, 1)
    if (debug) { print(paste("cur1p is", cur.prob.p)) }
  	
    for (pred in 1:np) {
  	  if (u.p[pred] <= (1 - cur.prob.p[pred])) {
  	    y.pred[i, pred] <- thresh
  	  } else {
  	  	adj <- (u.p[pred] - (1 - cur.prob.p[pred])) / (cur.prob.p[pred])
  	  	if (debug) {cat("adj is", adj, "\n")}
  	  	y.pred[i, pred] <- qGP(p=adj, scale=cur.sig.p[pred], shape=cur.xi.p[pred])
  	  }

  	}
  	
  }

  #################################################:
  ##########      Plotting            #############:
  #################################################:
  
  if (i %% update == 0) {
    if (debug) {cat("mean y is ", mean(y) )}
    if (iterplot) {
      par(mfrow = c(3, 2))
      
      plot(keep.beta.prob[1:i, 1], type="l")
      if (p1 > 1) {
      	plot(colMeans(keep.beta.prob[1:i, ]))
      } else {
        boxplot(keep.beta.prob[(i-999):i])
      }
      
      plot(exp(keep.beta.sig[1:i, 1]), type="l")
      if (p2 > 1) {
        plot(colMeans(exp(keep.beta.sig[1:i, ])))
      } else {
        boxplot(exp(keep.beta.sig[(i-999):i]))
      }
      
      plot(keep.beta.xi[1:i, 1], type="l")
      if (p3 > 1) {
        plot(colMeans(keep.beta.xi[1:i, ]))
      } else {
        boxplot(keep.beta.xi[(i-999):i])
      }

    }
    
  	cat("    Iter", i, "completed \n")
  }     
  
  
  }  # end iters 

  results <- list(beta.prob=keep.beta.prob, 
                  beta.sig=keep.beta.sig, 
                  beta.xi=keep.beta.xi, 
                  #dev=dev,
                  yp = y.pred)
                  
  return(results)
}