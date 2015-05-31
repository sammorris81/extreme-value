rm(list=ls())
load("simdata.RData")
ns <- dim(y)[1]
nt <- dim(y)[2]
nsettings <- dim(y)[4]
obs <- c(rep(T, 100), rep(F, 44))



setting.title <- c("Gaussian", "Symmetric-t (K = 1)",
                   "Symmetric-t (K = 5)",
                   paste("Skew-t (K = 1, ", expression(lambda), "= 3)"),
                   paste("Skew-t (K = 5, ", expression(lambda), "= 3)"),
                   "Max-stable", "Transform below T")
methods <- c("Gaussian", "Skew-t, K = 1, T = q(0.0)", 
             "Sym-t, K = 1, T = q(0.8)",
             "Skew-t, K = 5, T = q(0.0)", "Sym-t, K = 5, T = q(0.8)", 
             "Max-stable")
prefix <- c("gaus", "symt1", "symt5", "skewt1", "skewt5", "maxstab", "trans")
methodplot <- c(1, 6, 2, 3, 4, 5)

# reconstructing the time series

for (setting in 1:7) {
  for (method in methodplot) {
    set <- 1
    resultsfile <- paste("./results/",setting, "-", method, "-", set, ".RData", 
                         sep="")
    load(file=resultsfile)
    y.val <- y[!obs, , set, setting]
    s.val <- s[!obs, ]
    site <- 1
    if (method == 6) {
      yp <- fit.1[[1]]$yp[10001:20000, , site]  # method 6 transposes yp
    } else {
      yp <- fit.1$yp[, site, ]
    }
    pred.med <- apply(yp, 2, quantile, probs=0.5)
    pred.upp <- apply(yp, 2, quantile, probs=0.975)
    pred.low <- apply(yp, 2, quantile, probs=0.025)
    
    cover <- which(y.val[site, ] < pred.upp & y.val[site, ] > pred.low)
    ylim <- range(pred.med, pred.upp, pred.low, y.val[site, ])
    
    # plot the time series
    xplot <- seq(1 ,50)
    plot(pred.med, type="l", ylim=ylim, 
         main=paste("Data:", setting.title[setting], "\n Method:", 
                    methods[method]),
         xlab=paste("Time series at site", site), ylab="")
    lines(pred.upp, lty=2)
    lines(pred.low, lty=2)
    points(xplot[cover], y.val[site, cover], pch=21, 
           bg="dodgerblue1", col="dodgerblue4")
    points(xplot[-cover], y.val[site, -cover], pch=21, 
           bg="firebrick1", col="firebrick4")
    
    dev.print(device = pdf, paste(prefix[setting], "-", method, ".pdf", sep=""))
  }
}

prob.exceed <- function(mcmcoutput, thresh) {
  return (mean(mcmcoutput > thresh))
}

# visualizing the Brier score
for (setting in 1:7) {
  for (method in methodplot) {
    set <- 1
    resultsfile <- paste("./results/",setting, "-", method, "-", set, ".RData", 
                         sep="")
    load(file=resultsfile)
    y.val <- y[!obs, , set, setting]
    s.val <- s[!obs, ]
    thresh <- quantile(y.val, probs=0.95, na.rm=T)
    site   <- order(apply(y.val > thresh, 1, sum), decreasing = TRUE)[1]
    if (method == 6) {
      yp <- fit.1[[1]]$yp[10001:20000, , site]  # method 6 transposes yp
    } else {
      yp <- fit.1$yp[, site, ]
    }
    
    prob <- apply(yp, 2, prob.exceed, thresh=thresh)
    above <- which(y.val[site, ] > thresh)
    below <- which(y.val[site, ] <= thresh)
    
    # plot the probabilities
    xplot <- seq(1 ,50)
    plot(prob, type="l", ylim=c(0, 1), 
         main=paste("Data:", setting.title[setting], "\n Method:", 
                    methods[method]),
         xlab="Day", ylab="P[Y > q(0.95)]")
    points(above, rep(1, length(above)),
           pch=21, bg="dodgerblue1", col="dodgerblue4")
    points(below, rep(0, length(below)), 
           pch=21, bg="firebrick1", col="firebrick4")
    
    dev.print(device = pdf, paste("./plots/", prefix[setting], "-", method, "bs.pdf", sep=""))
  }
}

# Brier scores for q(0.98)
for (setting in 1:7) {
  for (method in methodplot) {
    set <- 1
    resultsfile <- paste("./results/",setting, "-", method, "-", set, ".RData", 
                         sep="")
    load(file=resultsfile)
    y.val <- y[!obs, , set, setting]
    s.val <- s[!obs, ]
    thresh <- quantile(y.val, probs=0.98, na.rm=T)
    site   <- order(apply(y.val > thresh, 1, sum), decreasing = TRUE)[1]
    if (method == 6) {
      yp <- fit.1[[1]]$yp[10001:20000, , site]  # method 6 transposes yp
    } else {
      yp <- fit.1$yp[, site, ]
    }
    
    prob <- apply(yp, 2, prob.exceed, thresh=thresh)
    above <- which(y.val[site, ] > thresh)
    below <- which(y.val[site, ] <= thresh)
    
    # plot the probabilities
    xplot <- seq(1 ,50)
    plot(prob, type="l", ylim=c(0, 1), 
         main=paste("Data:", setting.title[setting], "\n Method:", 
                    methods[method]),
         xlab="Day", ylab="P[Y > q(0.98)]")
    points(above, rep(1, length(above)),
           pch=21, bg="dodgerblue1", col="dodgerblue4")
    points(below, rep(0, length(below)), 
           pch=21, bg="firebrick1", col="firebrick4")
    
    dev.print(device = pdf, paste("./plots/", prefix[setting], "-", method, "bs.pdf", sep=""))
  }
}