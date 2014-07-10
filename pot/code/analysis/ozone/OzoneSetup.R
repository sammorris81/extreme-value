library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

#### Setup from Brian
load("OzoneData.RData")
y <- CMAQ_OZONE$y
x <- CMAQ_OZONE$x
s <- CMAQ_OZONE$s
l <- CMAQ_OZONE$poly

#Bounds of the subset for initial analysis
xb<-c(1201,1850)
yb<-c(-450,-250)
AL <- 1:26
MD1 <- 27:33
FL <- 37:45
GA <- 46:69
IL <- 70:76
IN <- 77:104
KY <- 105:135
MD2 <- 136:150
MS <- 151:155
NC <- 155:196
OH <- 197:223
PA <- 224:231
SC <- 232:252
TN <- 253:274

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

beta.init <- rep(0, dim(X)[3])
tau.init <- 1

save.image(file="cv-setup.RData")

# setting 1: gaussian
# setting 2: t0
# setting 3: t5
# setting 4: t0 T = 0.9
# setting 5: t5 T = 0.9

# chi plot
# bin information
d <- as.vector(dist(s))
j <- 1:54
i <- 2:55
ij <- expand.grid(i, j)
ij <- ij[(ij[1] > ij[2]), ]
sites <- cbind(d, ij)
dist <- rdist(s)
diag(dist) <- 0

bins <- seq(0, 600, 50)
bins <- c(bins, 900)

probs <- c(0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
threshs <- quantile(y, probs=probs, na.rm=T)
exceed <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(y[, t] > threshs[thresh])
    n.na <- sum(is.na(y[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(y[inbin, t]))
        acc[b] <- acc[b] + sum(y[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed[, thresh] <- exceed.thresh
}

par(mfrow=c(1, 2))
xplot <- (0:12) + 0.5
plot(xplot, exceed[, 1], type="b", ylim=c(0, 0.75), ylab="exceed", xaxt="n", xlab="bin distance", pch=1, lty=3, main="chi-plot ozone")
axis(1, at=0:12, labels=bins[-14])
for (line in 2:9) { 
	lines(xplot, exceed[, line], lty=3)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")

# Making the chi plot
# Residuals after lm
y.lm <- y[, 1]
x.lm <- X[, 1, ]
for (t in 2:nt) {
  y.lm <- c(y.lm, y[, t])
  x.lm <- rbind(x.lm, X[, t, ]) 
}
ozone.lm <- lm(y.lm ~ x.lm)
ozone.res <- residuals(ozone.lm)
ozone.int <- ozone.lm$coefficients[1]
ozone.beta1 <- ozone.lm$coefficients[3]
ozone.beta2 <- ozone.lm$coefficients[4]
ozone.cmaq <- ozone.lm$coefficients[5]
ozone.beta <- c(ozone.int, ozone.beta1, ozone.beta2, ozone.cmaq)
res <- matrix(NA, ns, nt)
for (t in 1:nt) {
  res[, t] <- y[, t] - X[, t, ] %*% ozone.beta
}

d <- as.vector(dist(s))
j <- 1:54
i <- 2:55
ij <- expand.grid(i, j)
ij <- ij[(ij[1] > ij[2]), ]
sites <- cbind(d, ij)
dist <- rdist(s)
diag(dist) <- 0

bins <- seq(0, 600, 50)
bins <- c(bins, 900)

probs <- c(0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
threshs <- quantile(res, probs=probs, na.rm=T)
exceed.res <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))

for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(res[, t] > threshs[thresh])
    n.na <- sum(is.na(res[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(res[inbin, t]))
        acc[b] <- acc[b] + sum(res[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed.res[, thresh] <- exceed.thresh
}

xplot <- (0:12) + 0.5
plot(xplot, exceed.res[, 1], type="b", ylim=c(0, 0.75), ylab="exceed", xaxt="n", xlab="bin distance", pch=1, lty=3, main="chi-plot ozone residuals")
axis(1, at=0:12, labels=bins[-14])
for (line in 2:9) { 
	lines(xplot, exceed.res[, line], lty=3)
	points(xplot, exceed.res[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")
