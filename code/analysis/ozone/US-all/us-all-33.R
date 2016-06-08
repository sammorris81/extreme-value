source("./package_load.R", chdir = TRUE)

setting <- 33
method <- "t"
nknots <- 6
keep.knots <- FALSE
threshold <- 0
tau.init <- 0.05
thresh.quant <- FALSE
skew <- TRUE
outputfile <- paste("results/us-all-", setting, ".RData", sep="")

start <- proc.time()

fit <- vector(mode="list", length=2)

for(val in 1:2){

	set.seed(setting*100 + val)

	cat("CV", val, "started \n")
	val.idx <- cv.lst[[val]]
	y.o <- Y[-val.idx,]
	X.o <- X[-val.idx, , ]
	S.o <- S[-val.idx,]

	y.p <- Y[val.idx,]
	X.p <- X[val.idx, , ]
	S.p <- S[val.idx,]

	tic.set <- proc.time()
	fit[[val]] <- mcmc(y=y.o, s=S.o, x=X.o, x.pred=X.p, s.pred=S.p,
                     temporalw=FALSE, temporalz=FALSE, temporaltau=FALSE,
                     method=method, skew=skew, keep.knots=keep.knots,
                     thresh.all=threshold, thresh.quant=thresh.quant, nknots=nknots,
                     iters=30000, burn=25000, update=500, iterplot=FALSE,
                     beta.init=beta.init, tau.init=tau.init, gamma.init=0.5,
                     rho.init=1, rho.upper=5, nu.init=0.5, nu.upper=10,
                     min.s=c(-2.25, -1.55), max.s=c(2.35, 1.30))
	toc.set <- proc.time()
	time.set <- (toc.set - tic.set)[3]

	elap.time.val <- (proc.time() - start)[3]
	avg.time.val <- elap.time.val / val
	cat("CV", val, "finished", round(avg.time.val, 2), "per dataset \n")
	save(fit, file=outputfile)
}
