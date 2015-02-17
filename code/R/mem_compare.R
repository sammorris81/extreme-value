library(inline)
library(SpatialTools)
library(fields)
library(microbenchmark)
# takes a mean and precision matrix and return a vector of conditional means and sds for specified rows
g.R <- function(d) {
  g <- apply(d, 1, which.min)
  results <- list(g=g)
  return(results)
}

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

s.temp <- cbind(runif(380, 0, 10), runif(380, 0, 10))
knots.temp <- cbind(runif(5, 0, 10), runif(5, 0, 10))
d.temp <- rdist(s.temp, knots.temp)
g.R(d=d.temp)
g.Rcpp(d=d.temp)

microbenchmark(g.R(d=d.temp), 
               g.mem(d=d.temp), times=1000)