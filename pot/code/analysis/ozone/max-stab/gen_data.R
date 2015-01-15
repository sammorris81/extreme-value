#############################################################:
###       FUNCTIONS TO GENERATE DATA FROM THE MODEL       ###:
#############################################################:


samp<-function(ns,ny,nk,x,beta,dw2,thresh=90){

    B<-array(0,c(ny,ns,5))
    for(j in 1:5){
      B[,,j]<-make.B(x,beta[,j],j)
    }
    alpha   <-B[1,1,4]
    gamma   <-B[1,1,5]

    tau<-Y<-matrix(0,ny,ns)
    for(t in 1:ny){
        a       <- rps(nk,alpha)
        A       <- a2AP(dw2,a,alpha,gamma)
        U       <- rGEV(ns,1,alpha,alpha)
        X       <- A*U
        tau[t,] <- exp(-1/X)
        Y[t,]   <- qGPD(tau[t,],B[t,,],thresh)
    }
Y}

rps<-function(n,alpha){
  #### PS(alpha) generation as given by Stephenson(2003)
  unif <- runif(n)*pi
  stdexp.ps <- rexp(n,1)
  logs <- (1-alpha)/alpha * log(sin((1-alpha)*unif)) + 
          log(sin(alpha*unif)) - (1-alpha)/alpha * log(stdexp.ps) - 
          1/alpha * log(sin(unif))
return(exp(logs))}

a2AP<-function(dw2,a,alpha,gamma){
    rho2 <- exp(gamma)^2
    fac  <- exp(-0.5*dw2/rho2)
    FAC  <- sweep(fac,1,rowSums(fac),"/")
    W    <- FAC^(1/alpha)
return((W%*%a)^alpha)}  


qGPD<-function(tau,B,thresh){

  if(is.matrix(B)){
    prob<-B[,1]
    sig<-B[,2]
    xi<-B[,3]
  }

  if(!is.matrix(B)){
    prob<-B[1]
    sig<-B[2]
    xi<-B[3]
  }
  
  tau2<-1-(tau-prob)/(1-prob)
  Y<-thresh+ifelse(tau<prob,0,sig*(tau2^(-xi)-1)/xi)
Y}


make.prob<-function(year,beta){
   eta<-beta[1]+year*beta[2]
   eta[eta>10]<-10
return(logit(eta))}

make.scale<-function(year,beta){
   eta<-beta[1]+year*beta[2]
   eta[eta>10]<-10
return(exp(eta))}

make.shape<-function(year,beta){
   eta<-beta[1]+year*beta[2]
return(eta)}

make.alpha<-function(year,beta){
   eta<-beta[1]+year*beta[2]
   eta[eta>10]<-10
return(logit(eta))}

make.gamma<-function(year,beta){
   eta<-beta[1]+year*beta[2]
   eta[eta>10]<-10
return(eta)}

make.B<-function(year,beta,type){
  B<-NA
  if(type==1){
    B<-make.prob(year,beta)
  }
  if(type==2){
    B<-make.scale(year,beta)
  }
  if(type==3){
    B<-make.shape(year,beta)
  }
  if(type==4){
    B<-make.alpha(year,beta)
  }
  if(type==5){
    B<-make.gamma(year,beta)
  }
return(B)} 

logit<-function(x){
  1/(1+exp(-x))
}

rGPD<-function(X,B,thresh){

  if(is.matrix(B)){
    prob<-B[,1]
    sig<-B[,2]
    xi<-B[,3]
  }

  if(!is.matrix(B)){
    prob<-B[1]
    sig<-B[2]
    xi<-B[3]
  }
  
  U<-exp(-1/X)
  U2<-1-(U-prob)/(1-prob)
  Y<-thresh+ifelse(U<prob,0,sig*(U2^(-xi)-1)/xi)
Y}

rGEV<-function(n,mu,sig,xi){
   tau<-runif(n)
   x<--1/log(tau)
   x<-x^(xi)-1
   x<-mu+sig*x/xi
x}
