In this zip file are all the necessary files to use Markov chain Monte Carlo 
(MCMC) methods to fit the space-time skew-t model for threshold exceedances. 

In the example.R file we include examples to generate the data using the 
rpotspatTS function in auxfunctions.R and fit the model using the mcmc function 
in mcmc.R. A description of the arguments for this function are included below.


#################################################################################
# MCMC
#
# Let ns be the number of sites, nt be the number of days, np be the number
# of parameters, and npred be the number of sites at which to make predictions.
# Finally, the model is written using sigma^2 as variance, but from a
# computation perspective, it's easier to store and use tau.
#
# y(ns, nt): data
# s(ns, 2): sites
# x(ns, nt, np): covariates
# s.pred(npred, 2): sites at which to make predictions
# x.pred(npred, nt, np): covariates for sites at which to make predictions
# min.s(2): smallest x and y
# max.s(2): largest x and y
# thresh.all(1): a scalar for value of the threshold (thresh.quant = FALSE) or
#                the quantile for the threshold (thresh.quant = TRUE)
# thresh.quant(1): boolean for whether thresholds given are quantiles (0, 1) or
#                  actual values from the data
# thresh.site.specific(1): boolean for site-specific threshold
# thresh.site(ns): vector of site-specific thresholds
# nknots(1): the number of knots to use in the MCMC
# fixknots(1): boolean on whether knot locations should be fixed
# keep.knots(1): should the MCMC return the knots in the output
# skew(1): boolean for whether a skew-model should be fit
# method(1): string ("t" or "gaussian") to indicate which process to fit
# cov.model(1): string ("matern" or "exponential")
# temporalw(1): boolean for TS on the knot location
# temporaltau(1): boolean for TS on the precision
# temporalz(1): boolean for TS on the skew term
#
# Starting values:
# beta.init(np): vector of starting values for beta
# tau.init(1 or nt): scalar or vector of starting values for precision
# tau.alpha.init(1): starting value of the a hyperparameter for precision
# tau.beta.init(1): starting value of the b hyperparameter for precision
# rho.init(1): starting value of the bandwidth for Matern Sigma
# nu.init(1): starting value of the smoothness for Matern Sigma
# gamma.init(1): starting value for strength of the nugget effect
# z.init(1 or nknots): starting values for z
# lambda.init(1): starting value for lambda
#
# Prior distributions:
# beta      ~iid N(beta.m, beta.s)
# tau.alpha ~ discrete(seq(tau.alpha.min, tau.alpha.max, by = tau.alpha.by))
# tau.beta  ~ gamma(tau.beta.a, tau.beta.b)
# rho       ~ uniform(0, rho.upper)
# lambda    ~ N(lambda.m, lambda.s)
#
# MCMC Settings:
# iters(1): number of iterations to run the MCMC
# burn(1): length of burnin period
# update(1): how often to update MCMC progress
# thin(1): how many observations to thin in the MCMC
# iterplot(1): boolean for whether iteration plots should be generated at an
#              update iteration for the MCMC
#
################################################################################