#########|#########|#########|#########|#########|#########|#########|#########
######################################################
########      THE MAIN MCMC FUNCTION       ###########
######################################################


Bayes_GPD <-function(
	y,								#data
	S,								#locations
	knots=NULL,						#knots
    X=NULL,							#covariates
    Sp=NULL, 						#Prediction locations
    Xp=NULL,						#Prediction covariates
    thresh=100,						#Threshold value
    init.beta=c(0,.1, -2.94),		#Initial GPD parameters - logscale, shape, pat
    init.alpha=.5,					#Initial alpha
    init.range=1,					#Initial range
    init.bw=1,						#Initial bandwidth
    init.logs=NULL,
    vary=c(F,F,F),					#Varying GPD parameters by location (logsig, xi, pat)
    pri.mn.range=-2,				#Range
    pri.sd.range=1,
    pri.mn.bw=0,					#Bandwidth
    pri.sd.bw=1,
    pri.var.a=0.1,					#a and b (gamma parameters for var)
    pri.var.b=0.1,
    pri.alpha.a=1,					#a and b (beta parameters for alpha)
    pri.alpha.b=1,
    #pri.mn.Phi.alpha=0,
    #pri.sd.Phi.alpha=0.25,
    pri.sd.beta=10,					#beta sd
    pri.mn.gev=c(0,0,-2.5),			    #mean for gev
    pri.sd.gev=c(1,0.25,0.5),		    #sd for gev params
    keep.samples=T,					#Retain samples
    nthreads = 2,
    iters=200,burn=50,update=50, nthin=1){

    start.time <- proc.time()

    if(is.null(knots)){
    	knots<-S
    }
    if(is.null(X)){
    	X<-cbind(1,S)
    }

	y[y <= thresh] <- thresh		#Everything below the threshold = threshold

    #### Some bookkeeping values
    n <- nrow(y)					#Number of sites
    nt <- ncol(y)					#Number of days
    nF <- nrow(knots)				#Number of knots
    p <- ncol(X)					#Number of covariates

    d<-rdist(S,S)					#Distance between locations
    diag(d)<-0
    dw2<-rdist(S,knots)^2			#Distance between locations and knots
    dw2[dw2<0.0001]<-0

    #### Initial values
    beta<-matrix(init.beta,n,3,byrow=T)	#GPD logscale, shape, and logit(pat) at each location
    mnb<-matrix(0,p,3)
    mnb[1,]<-colMeans(beta)			#Mean beta
    taub<-rep(1,3)					#Inv Variance for beta draws
    lograngeb<-log(init.range)		#Spatial range of GPD parameters
    logrange<-log(init.bw)			#Kernel bandwidth
    if(is.null(init.logs)){
    	logs<-matrix(2,nF,nt)		#log(A): A is basis function coefs
    } else {
    	logs<-matrix(init.logs,nF,nt)
    }		
    u<-matrix(0.5,nF,nt)			#Used for random effect: X = U*theta
    alpha<-init.alpha				#alpha term for random effects
    pat<-plogis(beta[,3])			#proportion that is thresholded
    thresh.mat<-y <= thresh			#which data<thresh


    #### Quantities used in posterior
    Qb<-solve(exp(-d/exp(lograngeb)))	#Sigma inv
    fac<-make.kern(dw2,logrange)	#Non-standardized kernel function weights
    FAC<-stdKern(fac)				#Standard kernel function weights
    curll<-0			#Current loglikelihood
    theta<-make.theta(FAC,logs,alpha)	#Come up with random effects
    
    curll <- dgpdcondCPP(y, thresh, beta[,1], beta[,2], pat, logs, FAC, theta, 
      					 alpha, logd=T, random=F, nthreads, loop.index=0,
      					 logs.debug=logs, theta.debug=theta, keepers)
    
    #CREATE PLACES TO STORE THE OUTPUT
    keepers<-matrix(0,iters,7)
    colnames(keepers)<-c("alpha","spatial range of the GEV params",
                         "kernel bandwidth",
                         "prop threshold","GEV log scale site 1",
                         "GEV shape site 1","log likelihood")
    GEVs<-c("GEV log scale","GEV shape", "logit(pat)")
    SITEs<-paste("Site",1:n)
    beta.mn<-beta.var<-0*beta
    colnames(beta.mn)<-colnames(beta.var)<-GEVs
    rownames(beta.mn)<-rownames(beta.var)<-SITEs
    locs=NULL						#Used only when keeping samples
    if(keep.samples){
      locs<-array(0,c(iters,n,3))
      dimnames(locs)<-list(paste("sample",1:iters),SITEs,GEVs)
    }
    Yp<-Yp1<-Yp2<-locsp<-beta.mn.p<-beta.var.p<-NULL
    if(!is.null(Sp)){
       np<-nrow(Sp)
       SITEs<-paste("Pred site",1:np)
       TIMEs<-paste("Replication",1:nt)
       Yp1<-Yp2<-matrix(0,np,nt)
       colnames(Yp1)<-colnames(Yp2)<-TIMEs
       rownames(Yp1)<-rownames(Yp2)<-SITEs
       beta.mn.p<-beta.var.p<-matrix(0,np,3)   
       colnames(beta.mn.p)<-colnames(beta.var.p)<-GEVs
       rownames(beta.mn.p)<-rownames(beta.var.p)<-SITEs
       if(keep.samples){
         locsp<-array(0,c(iters,np,3))
         dimnames(locsp)<-list(paste("sample",1:iters),SITEs,GEVs)
         Yp<-array(0,c(iters,np,nt))
         dimnames(Yp)<-list(paste("sample",1:iters),SITEs,TIMEs)
         thetapred<-array(0,c(iters,np,nt))
         dimnames(thetapred)<-list(paste("sample", 1:iters),SITEs,TIMEs)
       }
       dw2p<-rdist(Sp,knots)^2
       dw2p[dw2p<0.0001]<-0
       if(sum(vary)>0){
        d12<-rdist(Sp,S)
        d22<-rdist(Sp,Sp)
        d12[d12<0.0001]<-0
        diag(d22)<-0
       }
    }
	# Place to store logs
	logs.samples <- matrix(NA, iters, nF) 

    #SETTINGS FOR THE M-H CANDIDATE DISTRIBUTION
    attb<-accb<-MHb<-c(0.1,0.02,0.1,0.01,0.02,0.02,1)
    MHu<-3
    atts<-accs<-MHs<-c(3,2,2,2,rep(1,10), rep(0.5,10))
    cuts<-seq(0,40,2)


    #START SAMPLING!
    for(i in 1:iters){
    	
      for(rep in 1:nthin){

	  ##########################################################
      ##############      Random effects S and U    ############
      ##########################################################
      level<-olds<-logs
	  
      for(l in 1:nF){
      	W <- FAC[,l]^(1/alpha)
      	level[l,] <- vapply(logs[l,], get.level, FUN.VALUE=1, cuts=cuts)
      	MH1 <- MHs[level[l,]]
      	canlogs <- rnorm(nt, logs[l,], MH1)
      	canlevel <- vapply(canlogs, get.level, FUN.VALUE=1, cuts)
      	MH2 <- MHs[canlevel]
      	
      	# Trick to update theta without full matrix multiplication
      	cantheta <- theta + cbind(W)%*%rbind((exp(canlogs) - exp(logs[l,])))
      	
      	# Sometimes the trick results in a negative theta due to numerical
      	# instability. If this is the case then we should just remake the theta 
      	# matrix. Should be a mostly negligible time addition.
      	# Could probably be optimized to only recalculate for a specific
      	# day and location. 
      	if(sum(cantheta <= 0) > 0){
      		canlogsfull <- logs
      		canlogsfull[l,] <- canlogs
      		cantheta <- make.theta(FAC,canlogsfull,canalpha)
      	}
      	
      	canll <- dgpdcondCPP(y, thresh, beta[,1], beta[,2], pat, canlogs, FAC, 
      						 cantheta, alpha, logd=T, random=F, nthreads, 
      						 loop.index=0, logs.debug=logs, 
      						 theta.debug=theta, mcmcpart="logs",
      						 keepers)
      						   
      	R <- apply((canll - curll), 2, sum)+
      			h1(canlogs, u[l,], alpha, log=T)-
      			h1(logs[l,], u[l,], alpha, log=T)+
      			dnorm(logs[l,], canlogs, MH2, log=T)-
      			dnorm(canlogs, logs[l,], MH1, log=T)
      	
     	mh.compare <- rexp(nt, 1)
      	logs[l,(-R<mh.compare)] <- canlogs[(-R<mh.compare)]
      	curll[,(-R<mh.compare)] <- canll[,(-R<mh.compare)]
      	theta[,(-R<mh.compare)] <- cantheta[,(-R<mh.compare)]
      	
      }
           		   	
	logs.samples[i,] <- logs[,1]
      for(j in 1:length(MHs)){
         accs[j]<-accs[j]+sum(olds[level==j]!=logs[level==j])
         atts[j]<-atts[j]+sum(level==j)
      }
	
	#### Tune the candidate for log(A)
      for(j in 1:length(atts)){
      	if(i<burn/2 & atts[j]>50){
        	if(accs[j]/atts[j]<0.3){MHs[j]<-MHs[j]*0.9}
        	if(accs[j]/atts[j]>0.6){MHs[j]<-MHs[j]*1.1}
           	accs[j]<-atts[j]<-0
      	}
      }
	
	  #### B from article
      canu<-rtnorm(u)
      R<-h1(logs,canu,alpha,log=T)-
         h1(logs,u,alpha,log=T)+
         dtnorm(u,canu)-
         dtnorm(canu,u)
        
      acc<-matrix(rexp((nt*nF), 1), nF, nt)
      u<-ifelse(-R < acc, canu, u)

      ##########################################################
      ##############             alpha              ############
      ##########################################################
      canalpha<-rnorm(1,alpha,MHb[4])    
      if(i>50 & canalpha>0 & canalpha<1){
         attb[4]<-attb[4]+1
         cantheta<-make.theta(FAC,logs,canalpha)
         canll<-curll
         canll <- dgpdcondCPP(y, thresh, beta[,1], beta[,2], pat, logs, FAC, 
         					  cantheta, canalpha, logd=T, random=F, nthreads, 
         					  loop.index=0,
         					  logs.debug=logs, theta.debug=theta,
         					  mcmcpart="alpha", keepers)  
         					  
         R<-sum(canll-curll)+
            sum(h1(logs,u,canalpha))-
            sum(h1(logs,u,alpha))+
            dbeta(canalpha,pri.alpha.a,pri.alpha.b,log=T)-
            dbeta(alpha,pri.alpha.a,pri.alpha.b,log=T)

         if(-R < rexp(1, 1)){
           		alpha<-canalpha;curll<-canll;theta<-cantheta;
           		accb[4]<-accb[4]+1
         }           
      } 
      
      ##########################################################
      ##############       KERNEL BANDWIDTH         ############
      ##########################################################

      attb[5]<-attb[5]+1
      canlogrange<-rnorm(1,logrange,MHb[5])
      canfac<-make.kern(dw2,canlogrange)
      canFAC<-stdKern(canfac)
      canll<-curll
      cantheta<-make.theta(canFAC,logs,alpha)
      
      canll <- dgpdcondCPP(y, thresh, beta[,1], beta[,2], pat, logs, canFAC, 
      					   cantheta, alpha, logd=T, random=F, nthreads, 
      					   loop.index=0,
      					   logs.debug=logs, theta.debug=theta,
      					   mcmcpart="bandwidth", keepers)
      					   
      R<-sum(canll-curll)+
         dnorm(canlogrange,pri.mn.bw,pri.sd.bw,log=T)-
         dnorm(logrange,pri.mn.bw,pri.sd.bw,log=T)
      if(-R < rexp(1, 1)){
          logrange<-canlogrange;fac<-canfac
          FAC<-canFAC;theta<-cantheta;curll<-canll
          accb[5]<-accb[5]+1
      }

      #########################################################
      #############          GEV PARAMETERS        ############
      #########################################################

      ## SPATIALLY VARYING PARAMETERS
      
      for(l in 1:3){
      	if(i>50 & vary[l]){
        	
        	Xb<-X%*%mnb[,l]

	        for(j in 1:n){
	          	VVV<-taub[l]*Qb[j,j]
	          	MMM<-taub[l]*Qb[j,j]*Xb[j]-
	               taub[l]*sum(Qb[-j,j]*(beta[-j,l]-Xb[-j]))
	
	            attb[l]<-attb[l]+1
	            canb<-beta[j,]
	            canb[l]<-rnorm(1,beta[j,l],MHb[l])
	            canpat <- plogis(canb[3])
                
                #### Only evaluate the likelihood if sigma and xi give us a valid
          		#### likelihood for the data.
                if(!(canb[2] < 0 & y > (thresh - exp(canb[1])/canb[2]))){
	                canll <- dgpdcondCPP(y[j,], thresh, canb[1], canb[2], canpat, 
	                					 logs, FAC, theta[j,], alpha, logd=T, 
	                					 random=F, nthreads, loop.index=1,
	                					 logs.debug=logs, theta.debug=theta,
	                					 mcmcpart="GEV", keepers)
		
		            R<-sum(canll-curll[j,])+
		               dnorm(canb[l],MMM/VVV,1/sqrt(VVV),log=T)-
		               dnorm(beta[j,l],MMM/VVV,1/sqrt(VVV),log=T)
		            if(-R < rexp(1, 1)){
		            		beta[j,]<-canb
		               		curll[j,]<-canll
		               		accb[l]<-accb[l]+1
		        	}
	            }
	    	}

        }
      }
      

      # SPATIALLY CONSTANT PARAMETERS
      for(l in 1:3){
      	if(i>50 & !vary[l]){
          attb[l]<-attb[l]+1
          canb<-beta
          canb[,l]<-rnorm(1,beta[1,l],MHb[l])         
          canb[,l]<-beta[1,l]+MHb[l]*rt(1,df=5)
          canll<-curll
          canpat <- plogis(canb[,3])
          
          #### Only evaluate the likelihood if sigma and xi give us a valid
          #### likelihood for the data.
          if(sum(canb[,2] < 0 & (y > (thresh - exp(canb[,1])/canb[,2]))) == 0){
	          canll <- dgpdcondCPP(y, thresh, canb[,1], canb[,2], canpat, logs, 
	          					   FAC, theta, alpha, logd=T, random=F, nthreads, 
	          					   loop.index=0,
	          					   logs.debug=logs, theta.debug=theta,
	          					   mcmcpart="GEV", keepers)					   
	               
	          R<-sum(canll-curll)+
	             dnorm(canb[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)-
	             dnorm(beta[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)
	          if(-R < rexp(1, 1)){
	             beta<-canb
	             curll<-canll
	             accb[l]<-accb[l]+1
	          }
	      }
      }}

      pat <- plogis(beta[,3])

      ##########################################################
      ##############    Spatial hyperparameters     ############
      ##########################################################

      if(sum(vary)>0){
        tXQ<-t(X)%*%Qb
        tXQX<-tXQ%*%X
      }

      for(l in 1:3){
      	if(vary[l]){
	        #MEAN 
	        VVV<-solve(taub[l]*tXQX+(1/pri.sd.beta^2)*diag(p))
	        MMM<-taub[l]*tXQ%*%beta[,l]
	        mnb[,l]<-VVV%*%MMM+t(chol(VVV))%*%rnorm(p)
	  
	        #VARIANCE
	        SS<-quad.form(Qb,beta[,l]-X%*%mnb[,l])
	        taub[l]<-rgamma(1,n/2+pri.var.a,SS/2+pri.var.b)
        }
      }

      #SPATIAL RANGE
      if(sum(vary)>0){
       attb[6]<-attb[6]+1
       canlograngeb<-rnorm(1,lograngeb,MHb[6])
       canQb<-solve(exp(-d/exp(canlograngeb)))
       R<-0.5*sum(vary)*logdet(canQb)-
          0.5*sum(vary)*logdet(Qb)+
          dnorm(canlograngeb,pri.mn.range,pri.sd.range,log=T)-
          dnorm(lograngeb,pri.mn.range,pri.sd.range,log=T)
       for(l in 1:3){
       	if(vary[l]){
           R<-R-0.5*taub[l]*quad.form(canQb,beta[,l]-X%*%mnb[,l])
           R<-R+0.5*taub[l]*quad.form(Qb,beta[,l]-X%*%mnb[,l])
        }
       }
       if(-R < rexp(1,1)){
          lograngeb<-canlograngeb;Qb<-canQb
          accb[6]<-accb[6]+1
       }
      }

      #########  TUNE THE CANDIDATE DISTRIBUTION  #######
      for(j in 1:length(accb)){
      	if(i<burn/2 & attb[j]>50){
        	if(accb[j]/attb[j]<0.3){MHb[j]<-MHb[j]*0.9}
        	if(accb[j]/attb[j]>0.6){MHb[j]<-MHb[j]*1.1}
         	accb[j]<-attb[j]<-0
        }
      }

     }#end nthin

	#KEEP TRACK OF STUFF:
     keepers[i,]<-c(alpha,exp(lograngeb),
                    exp(logrange),plogis(beta[1,3]),beta[1,1:2],sum(curll))
     if(i>burn){
        nnn<-iters-burn
        beta.mn<-beta.mn+beta/nnn
        beta.var<-beta.var+beta*beta/nnn
     }
     if(keep.samples){locs[i,,]<-beta}

     #MAKE PREDICTIONS AT NEW LOCATIONS
     if(!is.null(Sp)){
       YYY<-matrix(0,np,nt)
       facp<-make.kern(dw2p,logrange)
       FACp<-stdKern(facp)
       thetap<-make.theta(FACp,logs,alpha)^alpha

       bp<-matrix(beta[1,],np,3,byrow=T)
       for(j in 1:3){
       	if(vary[j]){
          RRR<-beta[,j]-X%*%mnb[,j]
          bp[,j]<-Xp%*%mnb[,j]+
                  proj.beta(RRR,d12,d22,Qb,taub[j],lograngeb)
        }
       }

	   patp <- plogis(bp[,3])
	  
       YYY <- rgpdcond(
       		nreps.f=nt, S.f=Sp, knots.f=knots, thresh.f=thresh, 
       		logsig.f=bp[,1], xi.f=bp[,2], alpha.f=alpha, logrange.f=logrange, 
       		pat.f=patp, mu.f=NULL, sd.f=NULL, theta.f=thetap)$y

       if(i>burn){
         Yp1<-Yp1+YYY/(iters-burn)
         Yp2<-Yp2+YYY*YYY/(iters-burn)
         beta.mn.p<-beta.mn.p+bp/(iters-burn)
         beta.var.p<-beta.mn.p+bp*bp/(iters-burn)
       }
       if(keep.samples){
         Yp[i,,]<-YYY
         locsp[i,,]<-bp
         thetapred[i,,]<-thetap
       }
     }

     #DISPLAY CURRENT VALUE:
     # if(i%%1000==0){
       # accr <- accb/attb
       # par(mfrow=c(3,2))
       # plot(keepers[1:i,1],ylab="alpha",xlab="iteration",type="l", 
       	# main=bquote("ACCR" == .(accr[4])))
       # plot(keepers[1:i,3],ylab="bandwidth",xlab="iteration",type="l",
        # main=bquote("ACCR" == .(accr[5])))
       # plot(keepers[1:i,4],ylab=bquote(paste((1 - pi), " site 1")),xlab="iteration",type="l",
        # main=bquote("ACCR" == .(accr[3])))
       # plot(keepers[1:i,5],ylab="GEV log scale site 1",xlab="iteration",type="l",
        # main=bquote("ACCR" == .(accr[1])))
       # plot(keepers[1:i,6],ylab="GEV shape site 1",xlab="iteration",type="l",
        # main=bquote("ACCR" == .(accr[2])))
       # plot(keepers[1:i,7],ylab="Log likelihood",xlab="iteration",type="l")
     # }

		if(i%%1000==0){
		# par(mfrow=c(3, 3))
		# for (l in 1:9){
			# plot((logs.samples[1:i,l]), type="l")
		# }
			elap.time <- (proc.time() - tic.hm)[3]
			avg.time <- elap.time/(i/1000)
	    	print(paste("    Iter", i, "finished: ", round(avg.time, 2), " per 1000 iters"))
		}
	
    }

    stop.time <- proc.time()
    
    #dev.print(device=pdf, file="iterationplot.pdf")

#Return output:

list(time=stop.time-start.time,
     beta.var=beta.var-beta.mn^2,
     beta.mn=beta.mn,
     beta.samples=locs,
     beta.var.pred=beta.var.p-beta.mn.p^2,
     beta.mn.pred=beta.mn.p,
     beta.samples.pred=locsp,
     Y.mn.pred=Yp1,
     Y.var.pred=Yp2-Yp1^2,
     Y.samples.pred=Yp,
     theta.samples.pred=thetapred,
     parameters=keepers,
	 loga=logs.samples)
}

#############################################################
########       INCLUSION EXCLUSION MATRIX         ###########
#############################################################

create.comb <- function(K){
	# Initialize the inclusion/exclusion matrix
	incl.excl.mat <- cbind(rep(0, K))
	coefs <- 1
	
	for(k in 1:K){
		comb.mat <- combn(K, k)
		incl.excl.temp <- matrix(0, K, choose(K, k))
		for(i in 1:choose(K, k)){
			incl.idx <- comb.mat[,i]
			incl.excl.temp[incl.idx, i] <- 1
		}
		add.sub <- (-1)^k
		coefs <- c(coefs, rep(add.sub, choose(K, k)))
		incl.excl.mat <- cbind(incl.excl.mat, incl.excl.temp)
	}
	incl.excl.list <- list(coefs=coefs, incl.excl.mat=incl.excl.mat)
	return(incl.excl.list)
}

#############################################################
########   FUNCTIONS TO COMPUTE INITIAL VALUES    ###########
#############################################################

#### What sort of initial values are these?
get.inits.mu<-function(y,logsig=0,xi=.1){
  m<-median(y,na.rm=T)
  mu<-m-exp(logsig)*(log(2)^(-xi)-1)/xi
return(mu)}


######################################################
##########    FUNCTIONS USED FOR PREDICTION  #########
######################################################

proj.beta<-function(B,d12,d22,S11inv,tau,logrho){
   #B<-n-vector of observed beta (minus the mean)
   ns<-nrow(d22)
   rho<-exp(logrho)

   S22<-exp(-d22/rho)/tau
   S12<-exp(-d12/rho)/tau
   S11inv<-S11inv*tau

   P2<-S12%*%S11inv
   P1<-S22-S12%*%S11inv%*%t(S12)
   P1<-t(chol(P1))    

   Bnew<-P2%*%B+P1%*%rnorm(ns)
return(Bnew)}

######################################################
####  FUNCTION TO COMPUTE THE RANDOM EFFECTS  ########
######################################################

#### This function computes theta^(1/alpha). This makes it easier to iterate
#### over the theta loop when sampling for the random effects.
make.theta<-function(FAC.f,logs.f,alpha.f){
    #theta is n x t
    #alpha in (0,1)
    #logs = log(A) is nF x t
    #FAC is n x nF
    if(length(logs.f)==1){
       xxx<- FAC.f^(1/alpha.f) * logs.f}
    if(length(logs.f)>1){
       xxx<-(FAC.f^(1/alpha.f))%*%exp(logs.f)
	}
	return(xxx)
}  

#### Weights for the kernal basis functions
stdKern<-function(w,single=F){
  if(single){K<-w/sum(w)}   
  if(!single){K<-sweep(w,1,rowSums(w),"/")}
	return(K)
}  

#### Spatial kernel
make.kern<-function(d2,logrho){
   rho2<-exp(logrho)^2
   w<-exp(-0.5*d2/rho2)
	return(w)
}

######################################################
####         QUANTILE AND BRIER SCORES        ########
######################################################

#### Want to be able to calculate quantile and brier scores while the mcmc
#### program is running to save time and storage memory later.

#### To do: 
####    adjust arguments to function to bring in validation data
####    brier score function that returns a vector(validation thresholds)
####    quantile score function that returns a vector(quantiles)

######################################################
########            OTHER FUNCTIONS        ###########
######################################################

logd<-function(theta,v){
  sum(log(theta)-theta*v) 
}

#### Used to set the standard deviation for the candidate distribution
#### for the A terms in the random effect. When log(A) is large means
#### sd is smaller, and log(A) small means sd is larger. 
get.level<-function(logs,cuts){
    sum(logs>cuts)+1
}

logdet<-function(X){
  determinant(X)$modulus
}

trunc<-function(x,eps=.1){
  x<-ifelse(x<eps,eps,x)
  x<-ifelse(x>1-eps,1-eps,x)
x}

