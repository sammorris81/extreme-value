
##################################################:
######    Generate a data set and run MCMC   #####: 
##################################################:
library(fields)

#LOAD THE MCMC CODE:
source("...\\Bayes_GEV.R")  #You must change directories


#SPECIFY HOW TO GENERATE THE DATA
S<-expand.grid(1:7,1:7)           #DATA LOCATIONS
knots<-S                          #KNOT LOCATIONS
nt<-10                            #NUMBER OF YEARS
alpha<-.5                         #NUGGET EFFECT
mu<-S[,1]/10                      #GEV LOCATION
sig<-1                            #GEV SCALE
xi<-.1                            #GEV SHAPE
bw<-1                             #KERNEL BANDWIDTH


#GENERATE DATA FROM THE MODEL
y<-rgevspatial(nt,S,knots,mu,sig,xi,alpha,bw)


#ASSIGN COVARIATES IN THE PRIOR MEAN OF THE 
#SPATIALLY-VARYING GEV PARAMETERS
n<-nrow(S)
X<-as.matrix(cbind(1,S),n,3)


#WITHHOLD A TEST SET FOR PREDICTION
test<-runif(n)<.25
yp<-y[test,]
Sp<-S[test,]
Xp<-X[test,]


#FIT THE MODEL
iters<-5000   #This should be larger for a real analysis,
burn <-1000   #here its small for illustrative purposes.

fit<-Bayes_GEV(y=y[!test,],X=X[!test,],S=S[!test,],knots=knots,
               iters=iters,burn=burn,
               Xp=Xp,Sp=Sp,vary=c(T,F,F))


#PLOT THE POSTERIOR PREDICTIVE DISTRIBUTIONS FOR 5 YEARS
par(mfrow=c(3,2))
for(t in 1:5){
 boxplot(fit$Y.samples.pred[burn:iters,,t],ylim=range(yp),
         xlab="Spatial locations",ylab="Post pred dist")
 points(yp[,t],pch=19,col=2)
}

#Plot the true and estimated GEV location parameters:
plot(mu[!test],fit$beta.mn[,1],
     xlab="True GEV location",
     ylab="Post mean GEV location")
abline(0,1)


PROB<-yp
for(i in 1:nrow(yp)){for(j in 1:nt){
  PROB[i,j]<-mean(yp[i,j]>fit$Y.samples.pred[burn:iters,i,j])
}}
cover<-PROB>0.05 & PROB<0.95
print("Coverage of 90% prediction intervals")
print(mean(cover))
