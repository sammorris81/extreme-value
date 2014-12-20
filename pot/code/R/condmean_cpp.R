library(inline)
code <- '
  arma::vec m = Rcpp::as<arma::vec>(mn);
  arma::mat H = Rcpp::as<arma::mat>(prec);
  arma::vec tau = Rcpp::as<arma::vec>(taug);
  arma::vec r = Rcpp::as<arma::vec>(res);
  arma::vec inc = Rcpp::as<arma::vec>(include);
  int n = inc.size(); int idx;
  arma::rowvec H_temp;
  arma::colvec r_temp;
  NumericVector condmean(n);
  NumericVector condsd(n);
  double sd_i;
  for (int i = 0; i < n; i++) {
    idx = inc(i) - 1;
    H_temp = H.row(idx);
    r_temp = r % tau;
    H_temp.shed_col(idx);
    r_temp.shed_row(idx);
    sd_i = sqrt(1 / H(idx, idx)) / tau(i);
    condsd(i) = sd_i;
    condmean(i) = m(idx) - pow(sd_i, 2) * tau(i) * as_scalar(H_temp * r_temp);
  }
  return Rcpp::List::create(
    Rcpp::Named("cond.mn") = condmean, 
    Rcpp::Named("cond.sd") = condsd
  );
'
conditional.Rcpp <- cxxfunction(signature(mn="numeric", prec="numeric", taug="numeric", res="numeric", include="numeric"), code, plugin="RcppArmadillo")

conditional.mean <- function(mn, prec, res, taug=NULL, include=NULL) {
  if (is.null(taug)) {
    taug <- rep(1, length(mn))
  }
  if (is.null(include)) {
    include <- 1:length(mn)
  }
  
  results <- conditional.Rcpp(mn=mn, prec=prec, res=res, taug=taug, include=include)
}
