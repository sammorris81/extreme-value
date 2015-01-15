#Load the MCMC code

setwd("S:\\Documents\\SkewT\\")
source("gen_data.R")
source("S:\\Documents\\SkewT\\MCMC4MaxStable.R")
library(fields)


#Generate data

ny<-50  # of years
ns<-10  # of sites
nk<-5


#GPD parameters
thresh <- 90
prob0  <- c(-2,-.5)   # logit[P(below thresh)] = prob0[1] + prob0[2]*x
scale  <- c(1,.5)     # log(GPD scale) = prob0[1] + prob0[2]*x
shape  <- c(.1,0)     # GPD shape = shape[1] + shape[2]*x
alpha  <- c(0,0)      # logit(alpha) = alpha[1] + alpha[2]*x
gamma  <- c(0,0)      # log(kernel bandwidth) = gamma[1] + gamma[2]*x

B      <- cbind(prob0,scale,shape,alpha,gamma)

s      <- cbind(1:ns,0)
knots  <- cbind(seq(0,ns+1,length=nk),0)
x      <- (1:ny-ny/2)/10
x      <- matrix(x,ny,ns,byrow=FALSE)
dw2    <- rdist(s,knots)^2

set.seed(0820)
Y      <- samp(ns,ny,nk,x,B,dw2,thresh=thresh)


#Split into training and testing
test          <- rep(FALSE,ns)
test[c(1,ns)] <- TRUE
Yo            <- Y[,!test]
so            <- s[!test,]
xo            <- x[,!test]
Yp            <- Y[,test]
sp            <- s[test,]
xp            <- x[,test]
 

library(fields)
image.plot(1:ny,1:ns,Y,xlab="Year",ylab="Site")


#Fit the model

iters <- 500   # iterations
burn  <- 100   # burn-in

fit<-maxstable(Yo,xo,s=so,sp=sp,xp=xp,thresh=thresh,knots=knots,
               iters=iters,burn=burn,update=10,thin=10)


#Plot the results
X11()
boxplot(fit$samples[burn:iters,2,],
        main="Posterior dist of the time trends")
abline(h=0)

X11()
par(mfrow=c(2,1))
boxplot(fit$yp[,,1],outline=FALSE)
lines(Yp[,1],lwd=2,col=2)

boxplot(fit$yp[,,2],outline=FALSE)
lines(Yp[,2],lwd=2,col=2)
