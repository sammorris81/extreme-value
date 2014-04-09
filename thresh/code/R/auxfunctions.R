################################################################
# Arguments:
#   preds(yp, iters): mcmc predictions at validation
#                         locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
QuantScore <- function(preds, probs, validate){
  np <- length(validate)
  nprobs <- length(probs)
        
  # apply gives nprobs x n. looking to find each site's quantile over all
  # of the days.
  pred.quants <- apply(preds, 1, quantile, probs=probs, na.rm=T)
    
  scores.sites <- matrix(NA, nrow=nprobs, ncol=np)
    
  for (q in 1:nprobs) {
    diff <- pred.quants[q] - validate
    i <- ifelse(diff >= 0, 1, 0)
    scores.sites[q, ] <- 2 * (i - probs[q]) * diff
  }
    
  scores <- apply(scores.sites, 1, mean, na.rm=T)

  return(scores)
}

################################################################
# Arguments:
#   preds(yp, iters): mcmc predictions at validation
#                         locations
#   probs(nthreshs): sample quantiles for scoring
#   validate(np): validation data
#
# Returns:
#   list:
#     scores(nthreshs): a single brier score per threshold
################################################################
BrierScore <- function(preds, probs, validate){
  nthreshs <- length(probs)
  thresholds <- quantile(validate, probs=probs, na.rm=T)
    
  scores <- rep(NA, nthreshs)
  for (b in 1:nthreshs) {
    pat <- apply((preds > thresholds[b]), 1, mean)
    ind <- validate < thresholds[b]
    scores[b] <- mean((ind - pat)^2, na.rm=T)
  }
    
  return(scores)
}

#### From Brian - Not sure what it does yet
summarize <- function(y, digits=3) {
  p_ne_zero <- colMeans(y != 0)
  mean <- colMeans(y)
  sd <- apply(y, 2, sd)
  low90 <- apply(y, 2, quantile, 0.05)
  high90 <- apply(y, 2, quantile, 0.95)
   
  star <- ifelse(p_ne_zero > 0.5, 8, NA)
  star <- ifelse(p_ne_zero > 0.9, 88, star)
  star <- ifelse(p_ne_zero > 0.95, 888, star)


  X <- cbind(p_ne_zero, mean, sd, low90, high90)
  X <- round(X, digits)
  rownames(X) <- colnames(y)

X}

################################################################
# Arguments:
#   x(np): Vector of observations in (-Inf, Inf)
#
# Returns:
#   results(np): Inverse logit transformation of x
################################################################
expit <- function(x) {
  x <- as.vector(x)
  x <- ifelse(x > 10, 10, x)
  
  results <- exp(x) / (1 + exp(x))
  
  return(results)
}

################################################################
# Arguments:
#   x(np): Vector of observations in (0, 1)
#
# Returns:
#   results(np): Logit transformation of x
################################################################
logit <- function(x) {
  x <- as.vector(x)

  results <- log(x / (1 - x))
  
  return(results)
}

################################################################
# Arguments:
#   y(n): Vector of observations
#   scale(n): Scale parameters for the GPD
#   shape(n): Shape parameters for the GPD
#
# Returns:
#   llike(1): Sum of the current log likelihood
################################################################
dGP <- function(y, scale, shape){
  OK <- (y > 0) & ifelse (shape < 0, y < -scale / shape, TRUE)
  llike <- NA
  if (prod(OK) == 1) {
    llike <- -log(scale) -
             (1 + 1 / shape) * log(1 + shape * y / scale)
  }

  return(sum(llike))
}

################################################################
# Arguments:
#   p(n): vector of percentiles
#   scale(n): Scale parameters for the GPD
#   shape(n): Shape parameters for the GPD
#
# Returns:
#   llike(1): Sum of the current log likelihood
################################################################
qGP <- function(p, scale, shape){
  
  q <- scale * ((1 - p)^(-shape) - 1) / shape

  return(q)
}

