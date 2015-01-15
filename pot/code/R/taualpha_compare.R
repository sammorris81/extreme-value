library(inline)
library(SpatialTools)
library(fields)
library(microbenchmark)
# takes a mean and precision matrix and return a vector of conditional means and sds for specified rows
taualpha.R <- function(mmm, tau, tau.beta) {
  lll <- rep(0, length(mmm))
  for (l in 1:length(lll)) {
    lll[l] <- sum(dgamma(tau, mmm[l], tau.beta, log=T))
  }
  results <- list(lll=lll)
  return(results)
}

library(inline)
code <- '
  arma::vec m = Rcpp::as<arma::vec>(mmm);
  arma::mat t = Rcpp::as<arma::mat>(tau);
  double t_b = as<double>(tau_beta);
  int nm = m.size(); int nknots = t.n_rows; int nt = t.n_cols;
  NumericVector l(nm);
  for (int i = 0; i < nm; i++) {
    l(i) = 0;
    for (int j = 0; j < nknots; j++) {
      for (int k = 0; k < nt; k++) {
        l(i) += R::dgamma(t(j, k), m(i), 1/t_b, 1);
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("lll") = l
  );
'
taualpha.Rcpp <- cxxfunction(signature(mmm="numeric", tau="numeric", tau_beta="numeric"), code, plugin="RcppArmadillo")

taualpha.lll <- function(mmm, tau, tau.beta) {
  results <- taualpha.Rcpp(mmm=mmm, tau=tau, tau_beta=tau.beta)
  return(results)
}

lll <- mmm <- seq(0.1, 10, 0.1)
tau.temp <- matrix(rgamma(150, 3, 8), 5, 30)
tau.beta <- 8

taualpha.R(mmm=mmm, tau=tau.temp, tau.beta=tau.beta)
taualpha.Rcpp(mmm=mmm, tau=tau.temp, tau_beta=tau.beta)

microbenchmark(taualpha.R(mmm=mmm, tau=tau.temp, tau.beta=tau.beta), 
               taualpha.lll(mmm=mmm, tau=tau.temp, tau.beta=tau.beta), times=1000)