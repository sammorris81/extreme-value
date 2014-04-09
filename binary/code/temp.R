dw2.temp<-rdist(S.1,knots.1)^2			#Distance between locations and knots
dw2.temp[dw2.temp<0.0001]<-0
fac.temp<-make.kern(dw2.temp,logrange)	#Non-standardized kernel function weights
FAC.temp<-stdKern(fac.temp)		


y.1.temp <- y.1.1[,8,1]
K.temp <- sum(y.1.temp)
y.exceed.temp <- which(y.1.temp == 1)

n.temp <- length(y.1.temp)
incl.excl.temp <- matrix(1, n.temp, 2^K.temp)

all.comb <- create.comb(K.temp)

incl.excl.temp[y.exceed.temp,] <- all.comb$incl.excl.mat
coefs.temp <- all.comb$coefs

z.temp <- as.numeric( (1 + xi * X.hm.1 %*% beta.1)^(1/xi) )
inner.temp <- exp(1/alpha * (log(FAC.temp) - log(z.temp)))

pmf.temp <- rep(0, 2^K.temp)

for(k in 1:ncol(incl.excl.temp)){
	
	#### Selecting the proper locations to include in the inner sum
	inner.li.temp <- inner.temp * incl.excl.temp[,k]
	inner.l.temp <- colSums(inner.li.temp)^alpha
	
	#### This is not quite right yet. It shouldn't just be +/-
	#### it should be +/- in blocks
	pmf.temp[k] <- coefs.temp[k]*exp(-sum(inner.l.temp))
}

sum(pmf.temp)