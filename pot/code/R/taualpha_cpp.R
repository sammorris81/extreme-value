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