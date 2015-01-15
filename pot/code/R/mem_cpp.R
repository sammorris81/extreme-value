library(inline)
code <- '
  arma::mat d = Rcpp::as<arma::mat>(d_knots);
  int ns = d.n_rows; int nknots = d.n_cols;
  arma::uword index; double min;
  arma::rowvec d_temp(nknots);
  NumericVector g(ns);
  for (int i = 0; i < ns; i++) {
    d_temp = d.row(i);
    min = d_temp.min(index);
    g(i) = index + 1;
  }
  return Rcpp::List::create(
    Rcpp::Named("g") = g
  );
'
g.Rcpp <- cxxfunction(signature(d_knots="numeric"), code, plugin="RcppArmadillo")

g.mem <- function(d) {
  results <- g.Rcpp(d_knots=d)
}