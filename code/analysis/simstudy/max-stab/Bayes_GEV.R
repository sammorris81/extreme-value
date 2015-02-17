
######################################################
########      THE MAIN MCMC FUNCTION       ###########
######################################################


Bayes_GEV<-function(y,S,knots=NULL,
    X=NULL,
    Sp=NULL,Xp=NULL,
    init.beta=c(0,0,.0001),
    init.alpha=.25,
    init.range=1,
    init.bw=1,
    vary=c(T,F,F),
    pri.mn.range=-2,
    pri.sd.range=1,
    pri.mn.bw=0,
    pri.sd.bw=1,
    pri.var.a=0.1,
    pri.var.b=0.1,
    pri.alpha.a=1,
    pri.alpha.b=1,
    pri.sd.beta=10,
    pri.mn.gev=c(0,0,0),
    pri.sd.gev=c(10,1,0.25),
    keep.samples=T,
    iters=200,burn=50,update=50,nthin=1){

    start.time <- proc.time()
    library(fields)
    library(emulator)

    if(is.null(knots)){knots<-S}
    if(is.null(X)){X<-cbind(1,S)}


    #BOOKEEPING
    n<-nrow(y)
    nt<-ncol(y)
    nF<-nrow(knots)
    p<-ncol(X)

    d<-rdist(S,S)
    diag(d)<-0
    dw2<-rdist(S,knots)^2
    dw2[dw2<0.0001]<-0

    #INITIAL VALUES:

    beta<-matrix(init.beta,n,3,byrow=T)
    if(vary[1]){for(j in 1:n){
      beta[j,1]<-get.inits.mu(y[j,],beta[j,2],beta[j,3])
    }}
    mnb<-matrix(0,p,3)
    mnb[1,]<-colMeans(beta)
    taub<-rep(1,3)
    lograngeb<-log(init.range)
    logrange<-log(init.bw)
    logs<-matrix(2,nF,nt)
    u<-matrix(0.5,nF,nt)
    alpha<-init.alpha


    #COMPUTE QUANTITIES USED IN THE POSTERIOR
    Qb<-solve(exp(-d/exp(lograngeb)))
    fac<-make.kern(dw2,logrange)
    FAC<-stdKern(fac)
    curll<-matrix(0,n,nt)
    theta<-make.theta(FAC,logs,alpha)
    for(j in 1:n){
      curll[j,]<-loglike(y[j,],beta[j,1],beta[j,2],beta[j,3],theta[j,],alpha)
    }


    #CREATE PLACES TO STORE THE OUTPUT
    keepers<-matrix(0,iters,7)
    colnames(keepers)<-c("alpha","spatial range of the GEV params",
                         "kernel bandwidth",
                         "GEV loc site 1","GEV log scale site 1",
                         "GEV shape site 1","log likelihood")
    GEVs<-c("GEV location","GEV log scale","GEV shape")
    SITEs<-paste("Site",1:n)
    beta.mn<-beta.var<-0*beta
    colnames(beta.mn)<-colnames(beta.var)<-GEVs
    rownames(beta.mn)<-rownames(beta.var)<-SITEs
    locs=NULL
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


    #SETTINGS FOR THE M-H CANDIDATE DISTRIBUTION
    attb<-accb<-MHb<-c(0.1,0.02,0.02,0.01,0.02,0.02,1)
    MHu<-3
    atts<-accs<-MHs<-c(3,2,2,2,rep(1,10))
    cuts<-seq(0,15,2)


    #START SAMPLING!
    for(i in 1:iters){for(rep in 1:nthin){

      ##########################################################
      ##############      Random effects S and U    ############
      ##########################################################
      level<-olds<-logs

      for(t in 1:nt){ 
        v<-beta[,3]*exp(-beta[,2])*(y[,t]-beta[,1])+1
        v<-v^(-1/(alpha*beta[,3]))      
        ccc<-logd(theta[,t],v)
        for(l in 1:nF){
         W<-FAC[,l]^(1/alpha)
         level[l,t]<-get.level(logs[l,t],cuts)
         MH1<-MHs[level[l,t]]
         canlogs<-rnorm(1,logs[l,t],MH1)
         MH2<-MHs[get.level(canlogs,cuts)]

         cantheta<-theta[,t]+W*(exp(canlogs)-exp(logs[l,t]))
         canccc<-logd(cantheta,v)
         R<-canccc-ccc+
            h1(canlogs,u[l,t],alpha,log=T)-
            h1(logs[l,t],u[l,t],alpha,log=T)+
            dnorm(logs[l,t],canlogs,MH2,log=T)-
            dnorm(canlogs,logs[l,t],MH1,log=T)
         if(!is.na(exp(R))){if(runif(1)<exp(R)){
            logs[l,t]<-canlogs
            ccc<-canccc
            theta[,t]<-cantheta
         }}
        }
        curll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],theta[,t],alpha)
      }

      for(j in 1:length(MHs)){
         accs[j]<-accs[j]+sum(olds[level==j]!=logs[level==j])
         atts[j]<-atts[j]+sum(level==j)
      }
      for(j in 1:length(atts)){if(i<burn/2 & atts[j]>50){
        if(accs[j]/atts[j]<0.3){MHs[j]<-MHs[j]*0.9}
        if(accs[j]/atts[j]>0.6){MHs[j]<-MHs[j]*1.1}
           accs[j]<-atts[j]<-0
      }}

      canu<-rtnorm(u)
      R<-h1(logs,canu,alpha,log=T)-
         h1(logs,u,alpha,log=T)+
         dtnorm(u,canu)-
         dtnorm(canu,u)
      acc<-matrix(runif(nt*nF),nF,nt)
      u<-ifelse(acc<exp(R),canu,u)

      ##########################################################
      ##############              alpha             ############
      ##########################################################
      canalpha<-rnorm(1,alpha,MHb[4])    
      if(i>50 & canalpha>0 & canalpha<1){
         attb[4]<-attb[4]+1
         cantheta<-make.theta(FAC,logs,canalpha)
         canll<-curll  
         for(t in 1:nt){
           canll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta[,t],canalpha)
         }
         R<-sum(canll-curll)+
            sum(h1(logs,u,canalpha))-
            sum(h1(logs,u,alpha))+
            dbeta(canalpha,pri.alpha.a,pri.alpha.b,log=T)-
            dbeta(alpha,pri.alpha.a,pri.alpha.b,log=T)
         if(!is.na(exp(R))){if(runif(1)<exp(R)){
           alpha<-canalpha;curll<-canll;theta<-cantheta;
           accb[4]<-accb[4]+1
         }}           
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
      for(t in 1:nt){
        canll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta[,t],alpha)
      }
      R<-sum(canll-curll)+
         dnorm(canlogrange,pri.mn.bw,pri.sd.bw,log=T)-
         dnorm(logrange,pri.mn.bw,pri.sd.bw,log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
          logrange<-canlogrange;fac<-canfac
          FAC<-canFAC;theta<-cantheta;curll<-canll
          accb[5]<-accb[5]+1
      }}

      ##########################################################
      ##############          GEV PARAMETERS        ############
      ##########################################################

      ## SPATIALLY VARYING PARAMETERS
      for(l in 1:3){if(i>50 & vary[l]){
        Xb<-X%*%mnb[,l]

        for(j in 1:n){
          VVV<-taub[l]*Qb[j,j]
          MMM<-taub[l]*Qb[j,j]*Xb[j]-
               taub[l]*sum(Qb[-j,j]*(beta[-j,l]-Xb[-j]))

            attb[l]<-attb[l]+1
            canb<-beta[j,]
            canb[l]<-rnorm(1,beta[j,l],MHb[l])
            canll<-loglike(y[j,],canb[1],canb[2],canb[3],theta[j,],alpha)

            R<-sum(canll-curll[j,])+
               dnorm(canb[l],MMM/VVV,1/sqrt(VVV),log=T)-
               dnorm(beta[j,l],MMM/VVV,1/sqrt(VVV),log=T)
            if(!is.na(exp(R))){if(runif(1)<exp(R)){
               beta[j,]<-canb
               curll[j,]<-canll
               accb[l]<-accb[l]+1
            }}
          }

      }}

      ## SPATIALLY CONSTANT PARAMETERS
      for(l in 1:3){if(i>50 & !vary[l]){
          attb[l]<-attb[l]+1
          canb<-beta
          canb[,l]<-rnorm(1,beta[1,l],MHb[l])         
          canb[,l]<-beta[1,l]+MHb[l]*rt(1,df=5)
          canll<-curll
          for(t in 1:nt){
            canll[,t]<-loglike(y[,t],canb[,1],canb[,2],canb[,3],theta[,t],alpha)
          }
          R<-sum(canll-curll)+
             dnorm(canb[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)-
             dnorm(beta[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)
          if(!is.na(exp(R))){if(runif(1)<exp(R)){
             beta<-canb
             curll<-canll
             accb[l]<-accb[l]+1
          }}
      }}


      ##########################################################
      ##############    Spatial hyperparameters     ############
      ##########################################################

      if(sum(vary)>0){
        tXQ<-t(X)%*%Qb
        tXQX<-tXQ%*%X
      }

      for(l in 1:3){if(vary[l]){
        #MEAN 
        VVV<-solve(taub[l]*tXQX+(1/pri.sd.beta^2)*diag(p))
        MMM<-taub[l]*tXQ%*%beta[,l]
        mnb[,l]<-VVV%*%MMM+t(chol(VVV))%*%rnorm(p)
  
        #VARIANCE
        SS<-quad.form(Qb,beta[,l]-X%*%mnb[,l])
        taub[l]<-rgamma(1,n/2+pri.var.a,SS/2+pri.var.b)
      }}

      #SPATIAL RANGE
      if(sum(vary)>0){
       attb[6]<-attb[6]+1
       canlograngeb<-rnorm(1,lograngeb,MHb[6])
       canQb<-solve(exp(-d/exp(canlograngeb)))
       R<-0.5*sum(vary)*logdet(canQb)-
          0.5*sum(vary)*logdet(Qb)+
          dnorm(canlograngeb,pri.mn.range,pri.sd.range,log=T)-
          dnorm(lograngeb,pri.mn.range,pri.sd.range,log=T)
       for(l in 1:3){if(vary[l]){
           R<-R-0.5*taub[l]*quad.form(canQb,beta[,l]-X%*%mnb[,l])
           R<-R+0.5*taub[l]*quad.form(Qb,beta[,l]-X%*%mnb[,l])
       }}
       if(!is.na(exp(R))){if(runif(1)<exp(R)){
          lograngeb<-canlograngeb;Qb<-canQb
          accb[6]<-accb[6]+1
       }}
      }

      #########  TUNE THE CANDIDATE DISTRIBUTION  #######
      for(j in 1:length(accb)){if(i<burn/2 & attb[j]>50){
        if(accb[j]/attb[j]<0.3){MHb[j]<-MHb[j]*0.9}
        if(accb[j]/attb[j]>0.6){MHb[j]<-MHb[j]*1.1}
         accb[j]<-attb[j]<-0
      }}

     }#end nthin



     #KEEP TRACK OF STUFF:
     keepers[i,]<-c(alpha,exp(lograngeb),
                    exp(logrange),beta[1,],sum(curll))
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
       for(j in 1:3){if(vary[j]){
          RRR<-beta[,j]-X%*%mnb[,j]
          bp[,j]<-Xp%*%mnb[,j]+
                  proj.beta(RRR,d12,d22,Qb,taub[j],lograngeb)
       }}
       for(t in 1:nt){
          res<-rGEV(np,thetap[,t],alpha*thetap[,t],alpha)
          YYY[,t]<-bp[,1]+exp(bp[,2])*(res^(bp[,3])-1)/bp[,3]
       }

       if(i>burn){
         Yp1<-Yp1+YYY/(iters-burn)
         Yp2<-Yp2+YYY*YYY/(iters-burn)
         beta.mn.p<-beta.mn.p+bp/(iters-burn)
         beta.var.p<-beta.mn.p+bp*bp/(iters-burn)
       }
       if(keep.samples){
         Yp[i,,]<-YYY
         locsp[i,,]<-bp
       }
     }

     #DISPLAY CURRENT VALUE:
     if(i%%update==0){
       par(mfrow=c(3,2))
       plot(keepers[1:i,1],ylab="alpha",xlab="iteration",type="l")
       plot(keepers[1:i,3],ylab="bandwidth",xlab="iteration",type="l")
       plot(keepers[1:i,4],ylab="GEV loc site 1",xlab="iteration",type="l")
       plot(keepers[1:i,5],ylab="GEV log scale site 1",xlab="iteration",type="l")
       plot(keepers[1:i,6],ylab="GEV shape site 1",xlab="iteration",type="l")
       plot(keepers[1:i,7],ylab="Log likelihood",xlab="iteration",type="l")
     }

    }

    stop.time <- proc.time()

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
     parameters=keepers)}




#############################################################
########   FUNCTIONS TO COMPUTE INITIAL VALUES    ###########
#############################################################

get.inits.mu<-function(y,logsig=0,xi=.1){
  m<-median(y,na.rm=T)
  mu<-m-exp(logsig)*(log(2)^(-xi)-1)/xi
return(mu)}


######################################################
########   THE POSITIVE STABLE DENSITY    ###########
######################################################


h1<-function(logs,u,alpha,log=T){
   s<-exp(logs)
   psi<-pi*u
   c<-(sin(alpha*psi)/sin(psi))^(1/(1-alpha))
   c<-c*sin((1-alpha)*psi)/sin(alpha*psi)

   logd<-log(alpha)-log(1-alpha)-(1/(1-alpha))*logs+
         log(c)-c*(1/s^(alpha/(1-alpha)))+
         logs
   if(!log){logd<-exp(logd)}
logd}



######################################################
########           GEV FUNCTIONS           ###########
######################################################

loglike<-function(y,mu,logsig,xi,theta,alpha){
  missing<-is.na(y)
  theta<-theta^alpha
  mu.star<-mu+exp(logsig)*(theta^xi-1)/xi
  sig.star<-alpha*exp(logsig)*(theta^xi)
  xi.star<-alpha*xi
  ttt<-(1+xi.star*(y-mu.star)/sig.star)^(-1/xi.star)
  lll<- -log(sig.star)+(xi.star+1)*log(ttt)-ttt
  lll[missing]<-0
lll}

logd<-function(theta,v){
  sum(log(theta)-theta*v) 
}


rgevspatial<-function(nreps,S,knots,mu=1,sig=1,xi=1,alpha=.5,bw=1){
  library(evd)
    
  n<-nrow(S)
  nknots<-nrow(knots)

  d<-rdist(S,knots)
  d[d<0.0001]<-0
  w<-make.kern(d^2,log(bw))
  K<-stdKern(w)

  y<-matrix(0,n,nreps)
  for(t in 1:nreps){
    z<-K*t(rmvevd(n,alpha,d=nknots,mar=c(1,1,1)))
    X<-apply(z,1,max)
    y[,t]<-mu+sig*(X^xi-1)/xi  
  }  

return(y)}

######################################################
##########    FUNCTIONS USED FOR PREDICTION  #########
######################################################

rGEV<-function(n,mu,sig,xi){
   tau<-runif(n)
   x<--1/log(tau)
   x<-x^(xi)-1
   x<-mu+sig*x/xi
x}

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

make.theta<-function(FAC,logs,alpha){
    #theta is nxnF
    #s is nFxnt
    #alpha in (0,1)
    if(length(logs)==1){
       xxx<-(FAC^(1/alpha))*exp(logs)}
    if(length(logs)>1){
       xxx<-(FAC^(1/alpha))%*%exp(logs)}
xxx}  


stdKern<-function(w,single=F){
  if(single){K<-w/sum(w)}   
  if(!single){K<-sweep(w,1,rowSums(w),"/")}
K}  

make.kern<-function(d2,logrho){
   rho2<-exp(logrho)^2
   w<-exp(-0.5*d2/rho2)
w}



######################################################
########            OTHER FUNCTIONS        ###########
######################################################

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

rtnorm<-function(mn,sd=.25,fudge=0){
   upper<-pnorm(1-fudge,mn,sd)
   lower<-pnorm(fudge,mn,sd)
   if(is.matrix(mn)){
     U<-matrix(runif(prod(dim(mn)),lower,upper),dim(mn)[1],dim(mn)[2])
   }
   if(!is.matrix(mn)){
     U<-runif(length(mn),lower,upper)
   }
return(qnorm(U,mn,sd))}

dtnorm<-function(y,mn,sd=.25,fudge=0){
   upper<-pnorm(1-fudge,mn,sd)
   lower<-pnorm(fudge,mn,sd)
   l<-dnorm(y,mn,sd,log=T)-log(upper-lower)
return(l)}

