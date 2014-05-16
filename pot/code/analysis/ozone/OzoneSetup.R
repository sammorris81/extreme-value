library(fields)
library(geoR)
library(mvtnorm)

source('../simstudy/mcmc.R')
source('../simstudy/auxfunctions.R')

#### Setup from Brian
load("TempData.RData")
y <- CMAQ_OZONE$y
x <- CMAQ_OZONE$x
s <- CMAQ_OZONE$s
l <- CMAQ_OZONE$poly

#Bounds of the subset for initial analysis
xb<-c(1200,1850)
yb<-c(-500,-200)


#Plot the data

library(fields)
par(mfrow=c(2,2))

lim<-c(0,130)
plot(as.vector(x),as.vector(y),xlim=lim,ylim=lim,
     xlab="CMAQ output (ppb)",ylab="Monitor data (ppb)")
abline(0,1,lwd=2,col=2)

plot(s,axes=F,pch=19,main="Monitor locations",xlab="",ylab="")
lines(l)
abline(v=xb[1],col=2)
abline(v=xb[2],col=2)
abline(h=yb[1],col=2)
abline(h=yb[2],col=2)

image.plot(1:307,1:92,y,zlim=lim,xlab="Station",ylab="Day",main="Monitor data")

image.plot(1:307,1:92,x,zlim=lim,xlab="Station",ylab="Day",main="CMAQ data")


#Extract subset
NC<-(s[,1]>xb[1]) & (s[,1]<xb[2]) & (s[,2]>yb[1]) & (s[,2]<yb[2])
s<-s[NC,]
x<-x[NC,]
y<-y[NC,]
plot(s)
lines(l)

#########################################################################
#### Exclude if monitoring site is missing more than 50% of the 
#### observations
#########################################################################
excl <- which(rowMeans(is.na(y))>0.50)
s <- s[-excl,]
x <- x[-excl,]
y <- y[-excl,]

#### Rescale s so each dimension is in (0, 1)
s.scale <- matrix(NA, nrow=nrow(s), ncol=ncol(s))
s.scale[,1] <- (s[,1] - range(s[,1])[1])/(range(s[,1])[2] - range(s[,1])[1])
s.scale[,2] <- (s[,2] - range(s[,2])[1])/(range(s[,2])[2] - range(s[,2])[1])

#### Covariates including lat, long, and CMAQ
nt <- ncol(y)
ns <- nrow(y)
X <- array(1, dim=c(ns, nt, 4))
for(t in 1:nt){
	X[,t,2] <- s.scale[,1]	#Long
	X[,t,3] <- s.scale[,2]	#Lat
	X[,t,4] <- x[,t]
}

#### 5-fold cross validation
set.seed(2087)
cv.idx <- sample(1:nrow(s), nrow(s), replace=F)

cv.1 <- cv.idx[1:11]
cv.2 <- cv.idx[12:22]
cv.3 <- cv.idx[23:33]
cv.4 <- cv.idx[34:44]
cv.5 <- cv.idx[45:55]

cv.lst <- list(cv.1=cv.1, cv.2=cv.2, cv.3=cv.3, cv.4=cv.4, cv.5=cv.5)

ns <- nrow(y)
nt <- ncol(y)
iters <- 20000
burn <- 10000
probs <- c(0.80, 0.95, 0.99)
thresholds <- quantile(y, probs=probs, na.rm=T)

save.image(file="cv.RData")