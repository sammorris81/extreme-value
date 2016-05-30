fits <- list(fit.34, fit.39, fit.44, fit.60, fit.61, fit.62)
par(mfrow = c(3, 4))
intervals <- matrix(0, 6, 4)
for(i in 1:length(fits)) {
  intervals[i, 1:2] <- quantile(fits[[i]][[1]]$tau.beta, probs = c(0.025, 0.975))
  intervals[i, 3:4] <- quantile(fits[[i]][[2]]$tau.beta, probs = c(0.025, 0.975))
}

rm(list=ls())
setwd("/Volumes/sam-ext/skew-t/ozone-results")
done <- c(3:5, 7:9, 11:13, 15:17, 33:36, 38:41, 43:46, 51:74)
post.mean.a <- matrix(0, length(done), 2)
for (i in seq_along(done)) {
  setting <- done[i]
  filename <- paste("us-all-", setting, ".RData", sep = "")
  load(filename)
  post.mean.a[i, 1] <- mean(fit[[1]]$tau.alpha)
  post.mean.a[i, 2] <- mean(fit[[2]]$tau.alpha)
  print(paste("Setting ", setting, " complete", sep = ""))
}

mean(fit.34[[1]]$tau.alpha)
