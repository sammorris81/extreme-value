rm(list=ls())

load('simdata.RData')
load('debug.RData')
source('../../R/auxfunctions.R')
source('../../R/densities.R')

yactual <- full[,,15,2]

lll 		# All NaN

missing.obs	#none

thresh.obs 	#about 95%

mu.y 		# all reasonable

ns.v 		# 25 0 0 
r			# NaN on t=10, knot 1
r.mtx		# ok
res			# -Inf on day 10 site 2
str(Sig)	# PREC, logDET, SIG

PREC 	<- Sig$PREC		# Nothing unusual here
logDET 	<- Sig$logDET	# 0.4570, 0, 0
SIG		<- Sig$SIG		# Nothing unusual here


sig.r		# NaN
ss			# on day 10 NaN
thresh.mtx	# Nothing unusual
y.full		# -Inf
z			# Unusually large negative numbers

partition = matrix(1, 25, 10)

z[y.full==-Inf]
# -33.41145  -46.45216 -305.30466 -236.91665 -219.71632 -240.21376  -39.11015

z.impute

#### Reason why failing here:
# -Inf in y.full --> -Inf in res --> NaN in ss
thresh.mtx.z <- y2z(thresh.mtx, mu.y, r.mtx)
nt <- 10


for(i in 1:nrow(y.full)){
	
	
	thresh.days <- thresh.obs[i,]
	missing.days <- missing.obs[i,]
	           
	#### MVN impute
	S22 <- SIG[-i,]
	S22inv <- chol2inv(chol(S22[,-i]))
	Ez <- -SIG[i,-i] %*% S22inv %*% (z[-i,])
	Sz <- sqrt(SIG[i,i] - SIG[i,-i] %*% S22inv %*% SIG[-i,i])
	            
	z.impute <- rtnorm(mn=Ez, sd=Sz, lower=-Inf,
	                   upper=thresh.mtx.z[i,], fudge=0.000001)
	z.impute <- matrix(z.impute, 1, nt)
	z.missing <- rnorm(n=nt, mean=Ez, sd=Sz)
	z.missing <- matrix(z.missing, 1, nt)
	            
	y.impute <- z2y(z.impute, mu.y[i,], r.mtx[i,])
	y.missing <- z2y(z.missing, mu.y[i,], r.mtx[i,])
	            
	y.full[i, thresh.days] <- y.impute[thresh.days]
	y.full[i, missing.days] <- y.missing[missing.days]
}

SumSquares(res, Sig$PREC, partition, ns.v) # error in t(res.v[,t]) %*% PREC.v