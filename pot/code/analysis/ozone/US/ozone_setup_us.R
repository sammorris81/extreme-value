library(fields)
library(geoR)
library(mvtnorm)

rm(list=ls())
source('../../../R/mcmc.R')
source('../../../R/auxfunctions.R')

# Setup from Brian
load("ozone_data.RData")
S <- cbind(x[s[, 1]], y[s[, 2]])  # expands the grid of x, y locs where we have CMAQ 

# Exclude if site is missing more than 50% of its days
excl   <- which(rowMeans(is.na(Y)) > 0.50)
index  <- index[-excl]
Y      <- Y[-excl, ]
S      <- S[-excl, ]

image.plot(x, y, matrix(CMAQ[, 5], nx, ny), main="CMAQ output (ppb) - 2005-07-01 - no thin AQS sites")
points(S)       # Locations of monitoring stations
lines(borders)  # Add state lines

# only include sites from the eastern US
keep.these <- which(S[, 1] > 0)
index      <- index[keep.these]
Y          <- Y[keep.these, ]
S          <- S[keep.these, ]
cmaq       <- CMAQ[index, ]

# center and scale CMAQ data
cmaq <- (cmaq - mean(cmaq)) / sd(cmaq)

# add in intercept
nt <- ncol(Y)
X <- array(1, dim=c(nrow(cmaq), nt, 2))
for (t in 1:nt) {
  X[, t, 2] <- cmaq[, t]
}

# Plot CMAQ
image.plot(x, y, matrix(CMAQ[, 5], nx, ny), main="CMAQ output (ppb) - eastern US")
points(S)       # Locations of monitoring stations
lines(borders/1000)  # Add state lines

S <- S / 1000
x <- x / 1000
y <- y / 1000

#### 2-fold cross validation
cv.idx <- sample(nrow(S), nrow(S), replace=F)
cv.1 <- cv.idx[1:367]
cv.2 <- cv.idx[368:735]
cv.lst <- list(cv.1=cv.1, cv.2=cv.2)

beta.init <- 0
tau.init <- 1

save.image(file="cv-setup-us.RData")

# Remove non-US locations (for making maps - not done yet)
east <- s.p[, 1] > 2350
west <- s.p[, 1] < -2250
north <- s.p[, 2] > 1275
south <- s.p[, 2] < -1600

remove <- east | west | north | south


# setting 1: gaussian
# setting 2: t1
# setting 3: t5
# setting 4: t1 T = 0.9
# setting 5: t5 T = 0.9
# setting 6: skew-t1
# setting 7: skew-t1 T = 0.9
# setting 8: skew-t5

# setting 9 - 16: No CMAQ

# chi plot
# bin information
# d <- as.vector(dist(S))
# j <- 1:734
# i <- 2:735
# ij <- expand.grid(i, j)
# ij <- ij[(ij[1] > ij[2]), ]
# sites <- cbind(d, ij)

dist <- rdist(S)
diag(dist) <- 0

bins <- seq(0, 4, 0.5)

probs <- c(0.9, 0.95, 0.99)
threshs <- quantile(Y, probs=probs, na.rm=T)
exceed <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(Y[, t] > threshs[thresh])
    n.na <- sum(is.na(Y[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(Y[inbin, t]))
        acc[b] <- acc[b] + sum(Y[inbin, t] > threshs[thresh], na.rm=T)
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

# par(mfrow=c(1, 2))
xplot <- (1:(length(bins)-1)) - 0.5
plot(xplot, exceed[, 1], type="o", ylim=c(0, 0.5), ylab="exceed", xaxt="n", xlab="bin distance / 1000", pch=1, lty=3, main="chi-plot ozone")
axis(1, at=0:(length(bins)-2), labels=bins[-length(bins)])
for (line in 2:3) { 
	lines(xplot, exceed[, line], lty=line)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:3, pch=1:3, legend=probs, title="sample quants")

# Making the chi plot
# Residuals after lm
y.lm <- Y[, 1]
x.lm <- X[, 1]
for (t in 2:nt) {
  y.lm <- c(y.lm, y[, t])
  x.lm <- rbind(x.lm, X[, t, c(1, 2, 3)]) 
}
ozone.lm <- lm(y.lm ~ x.lm)
ozone.res <- residuals(ozone.lm)
ozone.int <- ozone.lm$coefficients[1]
ozone.cmaq <- ozone.lm$coefficients[5]
# ozone.beta <- c(ozone.int, ozone.beta1, ozone.beta2, ozone.cmaq)
ozone.beta <- c(ozone.int, ozone.beta1, ozone.beta2)
res <- matrix(NA, ns, nt)
for (t in 1:nt) {
  res[, t] <- y[, t] - X[, t, c(1, 2, 3)] %*% ozone.beta
}

d <- as.vector(dist(s))
j <- 1:85
i <- 2:86
ij <- expand.grid(i, j)
ij <- ij[(ij[1] > ij[2]), ]
sites <- cbind(d, ij)
dist <- rdist(s)
diag(dist) <- 0

bins <- seq(0, 600, 50)
bins <- c(bins, 900)

probs <- c(0.9, 0.95, 0.99)
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
plot(xplot, exceed.res[, 1], type="o", ylim=c(0, 0.75), ylab="exceed", xaxt="n", xlab="bin distance", pch=1, lty=1)
axis(1, at=0:12, labels=bins[-14])
for (line in 2:9) { 
	lines(xplot, exceed.res[, line], lty=line, pch=line, type="o")
}
legend("topright", lty=1:3, pch=1:3, legend=probs, title="sample quants")

# set.seed(2087)
# # nw: no obs > 75
# nw.these.l <- which((s[-exceed.75.these, 1] < 0) & (s[-exceed.75.these, 2] >= 0))
# n.nw.these.l <- length(nw.these.l) * 0.25
# nw.thin.l <- sample(nw.these.l, n.nw.these.l, replace=F)

# # nw: at least one obs > 75
# nw.these.h <- which((s[exceed.75.these, 1] < 0) & (s[exceed.75.these, 2] >= 0))
# n.nw.these.h <- length(nw.these.h) * 0.25
# nw.thin.h <- sample(nw.these.h, n.nw.these.h, replace=F)

# # sw: no obs > 75
# sw.these.l <- which((s[-exceed.75.these, 1] < 0) & (s[-exceed.75.these, 2] < 0))
# n.sw.these.l <- length(sw.these.l) * 0.25
# sw.thin.l <- sample(sw.these.l, n.sw.these.l, replace=F)

# # sw: at least one obs > 75
# sw.these.h <- which((s[exceed.75.these, 1] < 0) & (s[exceed.75.these, 2] < 0))
# n.sw.these.h <- length(sw.these.h) * 0.25
# sw.thin.h <- sample(sw.these.h, n.sw.these.h, replace=F)

# # ne: no obs > 75
# ne.these.l <- which((s[-exceed.75.these, 1] >= 0) & (s[-exceed.75.these, 2] >= 0))
# n.ne.these.l <- length(ne.these.l) * 0.25
# ne.thin.l <- sample(ne.these.l, n.ne.these.l, replace=F)

# # ne: at least one obs > 75
# ne.these.h <- which((s[exceed.75.these, 1] >= 0) & (s[exceed.75.these, 2] >= 0))
# n.ne.these.h <- length(ne.these.h) * 0.25
# ne.thin.h <- sample(ne.these.h, n.ne.these.h, replace=F)

# # se: no obs > 75
# se.these.l <- which((s[-exceed.75.these, 1] >= 0) & (s[-exceed.75.these, 2] < 0))
# n.se.these.l <- length(se.these.l) * 0.25
# se.thin.l <- sample(se.these.l, n.se.these.l, replace=F)

# # se: at least one obs > 75
# se.these.h <- which((s[exceed.75.these, 1] >= 0) & (s[exceed.75.these, 2] < 0))
# n.se.these.h <- length(se.these.h) * 0.25
# se.thin.h <- sample(se.these.h, n.se.these.h, replace=F)

# thin     <- c(nw.thin.l, nw.thin.h, sw.thin.l, sw.thin.h, ne.thin.l, ne.thin.h, se.thin.l, se.thin.h)

# index    <- index[thin]
# s        <- cmaq.s[index, ]  # only include the stratified sample of sites
# aqs.cmaq <- CMAQ[thin, ]  # CMAQ measurements at the AQS sites
# aqs      <- Y[thin, ]     # AQS measurements at the AQS sites


# Covariates including lat, long, and CMAQ
# ns <- nrow(aqs)
# nt <- ncol(aqs)
# X <- array(1, dim=c(ns, nt, 7))
# for(t in 1:nt){
  # X[, t, 2] <- s[, 1]    # Long
  # X[, t, 3] <- s[, 2]    # Lat
  # X[, t, 4] <- s[, 1]^2  # Long^2
  # X[, t, 5] <- s[, 2]^2  # Lat^2
  # X[, t, 6] <- s[, 1] * s[, 2]  # Interaction
  # X[, t, 7] <- aqs.cmaq[, t]   # CMAQ
# }

# X.scale <- array(1, dim=c(ns, nt, 7))
# for(t in 1:nt){
  # X[, t, 2] <- s.scale[, 1]    # Long
  # X[, t, 3] <- s.scale[, 2]    # Lat
  # X[, t, 4] <- s.scale[, 1]^2  # Long^2
  # X[, t, 5] <- s.scale[, 2]^2  # Lat^2
  # X[, t, 6] <- s.scale[, 1] * s.scale[, 2]  # Interaction
  # X[, t, 7] <- aqs.cmaq[, t]   # CMAQ
# }


#### 5-fold cross validation
# cv.idx <- sample(nrow(s.scale), nrow(s.scale), replace=F)

# cv.1 <- cv.idx[1:108]
# cv.2 <- cv.idx[108:216]
# cv.3 <- cv.idx[217:324]
# cv.4 <- cv.idx[325:432]
# cv.5 <- cv.idx[433:544]

# cv.lst <- list(cv.1=cv.1, cv.2=cv.2, cv.3=cv.3, cv.4=cv.4, cv.5=cv.5)

# Rescale s so each dimension is in (0, 1)
# S.scale      <- matrix(NA, nrow=nrow(S), ncol=ncol(S))
# S.scale[, 1] <- (S[, 1] - range(S[, 1])[1]) / (range(S[, 1])[2] - range(S[, 1])[1])
# S.scale[, 2] <- (S[, 2] - range(S[, 2])[1]) / (range(S[, 2])[2] - range(S[, 2])[1])

# # Statified sample of 50% to make covariance calculations easier
# # stratified by whether exceed 75
# # which sites have at no days that exceed 75, gives the row number of the index
# gthan.75.these <- which(rowMeans(aqs >= 75, na.rm=T) > 0)
# lthan.75.these <- which(rowMeans(aqs >= 75, na.rm=T) ==0)
# n.gthan.75 <- round(length(gthan.75.these) * 0.5)
# n.lthan.75 <- round(length(lthan.75.these) * 0.5)

# set.seed(2087)
# keep.gthan.these <- sample(gthan.75.these, n.gthan.75, replace=F)
# keep.lthan.these <- sample(lthan.75.these, n.lthan.75, replace=F)
# keep.these <- c(keep.gthan.these, keep.lthan.these)