library(inline)
library(SpatialTools)
library(fields)
library(microbenchmark)
# takes a mean and precision matrix and return a vector of conditional means and sds for specified rows
conditional.R <- function(mn, prec, res, include=NULL) {
  if (length(mn) != nrow(prec)) { print("dimension error")}
  if (is.null(include)) {
  	include = 1:length(mn)
  }
  cond.sd <- sqrt(1 / diag(prec))[include]
  cond.mean <- rep(NA, length(include))
  for (i in 1:length(include)) {
    idx <- include[i]
    cond.mean[i] <- mn[idx] - cond.sd[i]^2 * prec[idx, -idx] %*% res[-idx]
  }
  results <- list(cond.mn=cond.mean, cond.sd=cond.sd)
  return(results)
}

library(inline)
code <- '
  arma::vec m = Rcpp::as<arma::vec>(mn);
  arma::mat H = Rcpp::as<arma::mat>(prec);
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
    r_temp = r;
    H_temp.shed_col(idx);
    r_temp.shed_row(idx);
    sd_i = sqrt(1 / H(idx, idx));
    condsd(i) = sd_i;
    condmean(i) = m(idx) - pow(sd_i, 2) * as_scalar(H_temp * r_temp);
  }
  return Rcpp::List::create(
    Rcpp::Named("cond.mn") = condmean, 
    Rcpp::Named("cond.sd") = condsd
  );
'
conditional.Rcpp <- cxxfunction(signature(mn="numeric", prec="numeric", res="numeric", include="numeric"), code, plugin="RcppArmadillo")

conditional.mean <- function(mn, prec, res, include=NULL) {
  if (is.null(include)) {
    include = 1:length(mn)
  }
  
  results <- conditional.Rcpp(mn=mn, prec=prec, res=res, include=include)
}

s.temp <- cbind(runif(15, 0, 10), runif(15, 0, 10))
d.temp <- rdist(s.temp)
diag(d.temp) <- 0
cov <- simple.cov.sp(D=d.temp, sp.type="matern", sp.par=c(1, 1), error.var=0, smoothness=0.5, finescale.var=0)
prec <- solve(cov)
temp <- rnorm(15, 0, 1)
res <- rnorm(15, 0, 1) + temp
these <- which(temp < (max(temp) - 0.5))
conditional.R(mn = temp, prec=prec, res=res, include=these)
conditional.Rcpp(mn = temp, prec=prec, res=res, include=these)

microbenchmark(conditional.R(mn = temp, prec=prec, res=res, include=these), 
               conditional.mean(mn = temp, prec=prec, res=res, include=these), times=1000)