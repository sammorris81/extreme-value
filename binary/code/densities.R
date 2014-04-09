###############################################################################
###### pmf eval                                                          ######
###############################################################################

pspatbin <- function(alpha.f, FAC.f, xi.f, beta.f, X.f, y.f){
	
	#### A single day
	#### To do: Rewrite function to take all observations at once
	n <- length(y.f)
	K <- sum(y.f)
	
	#### Identify the indices that exceed the threshold
	y.exceed <- which(y.f == 1)
	
	#### Create inclusion/exclusion matrix
	#### Starting with a matrix that is only 1
	incl.excl <- matrix(1, n, 2^K)
	
	#### Now change entries for proper inclusion/exclusion
	incl.excl[y.exceed,] <- create.comb(K)	
	z <- (1 + xi.f * X.f %*% beta.f)^(1/xi.f)			#should be n x 1
	inner <- exp(1/alpha.f * ( log(FAC.f) - log(z) ) )	#should be n x L
	
	
}

######################################################
#### Positive Stable Density                      ####
#### From Reich and Shaby (2011)                  ####
######################################################

h1<-function(logs,u,alpha,log=T){
    s<-exp(logs)
    psi<-pi*u
    c<-(sin(alpha*psi)/sin(psi))^(1/(1-alpha))
    c<-c*sin((1-alpha)*psi)/sin(alpha*psi)
  
    logd<-log(alpha)-log(1-alpha)-(1/(1-alpha))*logs+
          log(c)-c*(1/s^(alpha/(1-alpha)))+logs
        
	if(!log){
		logd<-exp(logd)
  	}
  
  return(logd)
}

###############################################################################
###### cdf eval                                                          ######
###############################################################################



dtnorm<-function(y,mn,sd=.25,fudge=0){
	upper<-pnorm(1-fudge,mn,sd)
	lower<-pnorm(fudge,mn,sd)
	l<-dnorm(y,mn,sd,log=T)-log(upper-lower)

	return(l)
}


###############################################################################
###### generate r.v.s                                                    ######
###############################################################################

rgevcond <- function(nreps.f, S.f, knots.f, alpha.f=.5, logrange.f=1, theta.f=NULL){
		
	n<-nrow(S.f)
	nknots<-nrow(knots.f)
	
	d<-rdist(S.f,knots.f)
	d[d<0.0001]<-0
	w<-make.kern(d^2,logrange.f)
	
	#### K is n x nF
	K<-stdKern(w)
	
	y<-matrix(0,n,nreps.f)
	
	simsetup <- F
	
	if(is.null(theta.f)){
		simsetup <- T
		#### PS(alpha) generation as given by Stephenson(2003)
		unif <- runif(nknots*nreps.f, 0, 1) * pi
		stdexp <- rexp(nknots*nreps.f, 1)
		logs <- (1-alpha.f)/alpha.f * log(sin((1-alpha.f)*unif)) + 
				log(sin(alpha.f*unif)) - (1-alpha.f)/alpha.f * log(stdexp) - 
				1/alpha.f * log(sin(unif))
				
		#### logs is nF x nt
		logs <- matrix(logs,nknots,nreps.f)
		
		#### theta is n x nt
		theta.f <- make.theta(K, logs, alpha.f)
		theta.f <- theta.f^alpha.f
	}		
	
	for(t in 1:nreps.f){
		u <- rgev(n, 1, alpha.f, alpha.f)
		x <- u * theta.f[,t]
		y[,t] <- x
	}
		 
	if(simsetup){
		ret.f <- list(y=y, theta=theta.f, logs=logs)
	} else {
		ret.f <- list(y=y)
	}
	
	return(ret.f)
}

######################################################
#### Independent GEV                              ####
######################################################
rgev<-function(n,mu,sig,xi){
  tau<-runif(n)
  x<--1/log(tau)
  x<-x^(xi)-1
  x<-mu+sig*x/xi
  
  return(x)
}

######################################################
#### Truncated Normal                             ####
######################################################
rtnorm<-function(mn,sd=.25,fudge=0){
	upper<-pnorm(1-fudge,mn,sd)
	lower<-pnorm(fudge,mn,sd)
	if(is.matrix(mn)){
		U<-matrix(runif(prod(dim(mn)),lower,upper),dim(mn)[1],dim(mn)[2])
	}
	if(!is.matrix(mn)){
		U<-runif(length(mn),lower,upper)
	}
	return(qnorm(U,mn,sd))
}

######################################################
########            OTHER FUNCTIONS        ###########
######################################################

######################################################
#### GPD quantiles                             	  ####
######################################################

qgpdcond <- function(p.f, thresh.f, logsig.f, xi.f, theta.f, alpha.f, pat.f){
	
	nt <- ncol(theta.f)
	n <- nrow(theta.f)
	y <- matrix(NA,n,nt)
	
	if(length(xi.f) == 1){
		xi.f <- rep(xi.f, n)
	}
	
	if(length(logsig.f) == 1){
		logsig.f <- rep(logsig.f, n)
	}

	if(length(pat.f) == 1){
		pat.f <- rep(pat.f, n)
	}
	
	u <- qgev(p.f, 1, alpha.f, alpha.f)
	a <- matrix(NA, n, nt)
	sig <- exp(logsig.f)
	
	for(t in 1:nt){
		X <- theta.f[,t] * u
			
		#### Transform the data to Unif(0, 1)
	   	a[,t] <- exp(-1/X)
	   	#### Set threshold for all values below pi. Trying to avoid a loop 
	   	#### here. The idea is to set everything to its transformed value then
	   	#### reassign the y values that should be the threshold.
	   	
	   	#### Condition that we're above threshold. New percentile is based
	   	#### on a > pi
	   	a.cond <- (a[,t]-(1-pat.f))/pat.f
	   	valid <- a.cond > 0
	   	if(sum(valid)>0){
	   		y[valid,t] <- qgpd2(a.cond[valid], loc=thresh.f, sig[valid], xi.f[valid])
	   	}
	    	
	   	y[!valid,t] <- thresh.f
	}
	return(y)	
}

#### Rewriting this function to allow for a vector of shape parameters.
#### Original function required length(shape)=1, but this is not necessary
#### if you use ifelse statement.
qgpd2 <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, 
    lambda = 0) 
{
    if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
        1) 
        stop("`p' must contain probabilities in (0,1)")
    if (min(scale) < 0) 
        stop("invalid scale")
    if ((lambda < 0) || (lambda >= 1) || length(lambda) != 1) 
        stop("invalid lambda")
    if (any(p < lambda)) 
        stop("``p'' must satisfy ``p >= lambda''")
    if (lower.tail) 
        p <- 1 - p
    p <- p/(1 - lambda)
    ifelse(
    	shape == 0, 
        return(loc - scale * log(p)),
        return(loc + scale * (p^(-shape) - 1)/shape)
    )
}


#### Function to calculate the probability of being thresholded given theta
cond.dens.thresh <- function(alpha.f, pat.f, theta.f, log=F#,
		#logs.f, FAC.f, logrange.f, y.f, valid #### debug
	){
	py.f <- - theta.f * (-log(1-pat.f))^(1/alpha.f)
	
	if(!log){
		py.f <- exp(py.f)
	}
	return(py.f)
}


