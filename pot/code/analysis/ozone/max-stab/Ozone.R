#Load the MCMC code

load("../ozone_data.RData")
source("gen_data.R")
source("MCMC4MaxStable.R")
library(fields)
 
thresh<-75
Y <- t(Y)
Y <- ifelse(Y>thresh,Y,thresh)
X <- t(CMAQ[index,])
S <- cbind(x[s[,1]],y[s[,2]])

S <- S/1000
X <- (X-60)/15

#start with a subset
Y <- Y[,1:100]
X <- X[,1:100]
S <- S[1:100,]

#splint into training and testing
test <- rank(runif(100))<50
Yo   <- Y[,!test]
Xo   <- X[,!test]
So   <- S[!test,]
Yp   <- Y[,test]
Xp   <- X[,test]
Sp   <- S[test,]

knots<-So

iters <- 5000   # iterations
burn  <- 1000   # burn-in

fit<-maxstable(Yo,Xo,s=So,sp=Sp,xp=Xp,thresh=thresh,knots=knots,
               iters=iters,burn=burn,update=10,thin=1)

t<-2
boxplot(fit$yp[,t,],outline=FALSE)
lines(Yp[t,],lwd=2,col=2)





