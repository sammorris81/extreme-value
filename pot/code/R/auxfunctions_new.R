#########################################################################
# Arguments:
#   mn(nt): mean
#   sd(nt, nt): standard deviation
#   lower(1): lower truncation point (default=-Inf) 
#   upper(1): upper truncation point (default=Inf)
#	fudge(1): small number for numerical stability (used for lower bound)
#
# Returns:
#   y(nt): truncated normal data
#########################################################################
rTNorm <- function(mn, sd, lower=-Inf, upper=Inf, fudge=0){
  lower.u <- pnorm(lower, mn, sd)
  upper.u <- pnorm(upper, mn, sd)
  
  # replace <- ((mn / sd) > 5) & (lower == 0)
  # lower.u[replace] <- 0
  lower.u <- ifelse( mn / sd > 5 & lower == 0, 0, lower.u )
  U <- runif(length(mn), lower.u, upper.u)
  y <- qnorm(U, mn, sd)
  
  return(y)
}

CorFx <- function(d, alpha, rho, nu){

  library(geoR)    
  cor       <- alpha * matern(d, rho, nu)
  diag(cor) <- 1

return(cor)}


mem<-function(s,knots){
   library(fields)
   d<-rdist(s,knots)
   g<-apply(d,1,which.min)
g}