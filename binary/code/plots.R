###############################################################################
#### Basic setup                                                           ####
###############################################################################
#setwd('~/svn-repos/dissertation/extreme-value/code')
#setwd('S:/Documents/Dissertation/extreme-value/code')

library(fields)
library(emulator)
library(evd)
library(POT)
library(geoR)

load(file="./simdata.RData")

#### Model parameters                
alpha <- c(.3, .5, .7)  # Spatial correlation 
                            # (0 = spatial dep, 1 = spatial ind)     
logrange <- 0         				    # Kernel bandwidth 

for(set in 1:nsets){
	for(day in 1:nt){
		file1 <- paste("dataplots/sim1day",day,".pdf", sep="")
		file2 <- paste("dataplots/sim2day",day,".pdf", sep="")
		file3 <- paste("dataplots/sim3day",day,".pdf", sep="")
		dev.new(width=10, height=5)
		par(mfrow=c(1, 2))
		z <- matrix(y.1[,day,set], length(s2), length(s2))
		image.plot(s1, s2, z, 
			main=bquote(paste(
				alpha == .(alpha[1]), ", ", xi==.(xi), " I(Z > Thresh)")
			)
		)
		z <- matrix(latent.1[,day,set], length(s2), length(s2))
		image.plot(s1, s2, z, 
			main=bquote(paste(
				alpha == .(alpha[1]), ", ", xi==.(xi), " Latent Process")
			)
		)
		dev.print(file=file1, device=pdf)
		
		z <- matrix(y.2[,day,set], length(s2), length(s2))
		image.plot(s1, s2, z, 
			main=bquote(paste(
				alpha == .(alpha[2]), ", ", xi==.(xi), " I(Z > Thresh)")
			)
		)
		z <- matrix(latent.2[,day,set], length(s2), length(s2))
		image.plot(s1, s2, z, 
			main=bquote(paste(
				alpha == .(alpha[2]), ", ", xi==.(xi), " Latent Process")
			)
		)
		dev.print(file=file2, device=pdf)
		
		z <- matrix(y.3[,day,set], length(s2), length(s2))
		image.plot(s1, s2, z, 
			main=bquote(paste(
				alpha == .(alpha[3]), ", ", xi==.(xi), " I(Z > Thresh)")
			)
		)
		z <- matrix(latent.3[,day,set], length(s2), length(s2))
		image.plot(s1, s2, z, 
			main=bquote(paste(
				alpha == .(alpha[3]), ", ", xi==.(xi), " Latent Process")
			)
		)
		dev.print(file=file3, device=pdf)
		
		
		dev.off()
	}
}


#par(mfrow=c(2,2))

xplota <- seq(75, 100, by=0.01)
xplotb <- seq(100.01, 140, by=0.01)

yplot1a <- dnorm(xplota, mu[1], sd[1])
yplot1b <- pat[1] * dgpd(xplotb, thresh, exp(logsig[1]), xi[1])

plot(xplota, yplot1a, type="l", lty=1, xlim=c(75, 140))
lines(xplotb, yplot1b, type="l", lty=1)

yplot2a <- dnorm(xplota, mu[2], sd[2])
yplot2b <- pat[2] * dgpd(xplotb, thresh, exp(logsig[2]), xi[2])

lines(xplota, yplot2a, type="l", lty=2)
lines(xplotb, yplot2b, type="l", lty=2)

yplot3a <- dnorm(xplota, mu[3], sd[3])
yplot3b <- pat[3] * dgpd(xplotb, thresh, exp(logsig[3]), xi[3])

lines(xplota, yplot3a, type="l", lty=3)
lines(xplotb, yplot3b, type="l", lty=3)

yplot4a <- dnorm(xplota, mu[4], sd[4])
yplot4b <- pat[4] * dgpd(xplotb, thresh, exp(logsig[4]), xi[4])

lines(xplota, yplot4a, type="l", lty=4)
lines(xplotb, yplot4b, type="l", lty=4)

############################################################################

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