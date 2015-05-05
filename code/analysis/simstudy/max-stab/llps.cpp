// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat dPS_cpp_mat(arma::mat a, double alpha, arma::vec psi,
                      arma::vec bin_width, int threads) {

#ifdef _OPENMP
  if(threads > 0) {
    omp_set_num_threads(threads);
  }
#endif

  int ns = a.n_rows; int nt = a.n_cols;
  int nbins = psi.n_elem;
  int s; int t; int i;
  double integral; double llst; double logc; double logint;
  double ast; double psi_i;
  arma::mat ll(ns, nt); // arma::vec psi(nbins);

  // psi = dscal(nbins, PI, mid_points, 1);
  // psi = PI * mid_points;

#pragma omp parallel for \
  private(t, i, ast, llst, logc, logint, integral, psi_i) \
  shared(a, alpha, ll, psi) \
  schedule(static)
  for (s = 0; s < ns; s++) {
    for (t = 0; t < nt; t++) {
      ast = a(s, t);
      llst = log(alpha) - log(1 - alpha) - log(ast) / (1 - alpha);
      integral = 0;
      for (i = 0; i < nbins; i++) {
        psi_i = psi[i];
        logc = (log(sin(alpha * psi_i)) - log(sin(psi_i))) / (1 - alpha) +
          log(sin((1 - alpha) * psi_i)) - log(sin(alpha * psi_i));
        logint = logc - exp(logc) * pow(ast, (- alpha / (1 - alpha)));
        integral += exp(logint) * bin_width[i];
      }
      ll(s, t) = llst + log(integral);
    }
  }

  return ll;
}

// [[Rcpp::export]]
double dPS_cpp_sca(double ast, double alpha, arma::vec psi,
                   arma::vec bin_width, int threads) {

#ifdef _OPENMP
  if(threads > 0) {
    omp_set_num_threads(threads);
  }
#endif

  int nbins = psi.n_elem;
  int i;
  arma::vec integral(nbins);
  double llst; double psi_i; double logc; double logint;
  double ll;

  // psi = dscal(nbins, PI, mid_points, 1);

  llst = log(alpha) - log(1 - alpha) - log(ast) / (1 - alpha);
#pragma omp parallel for \
  private(psi_i, logc, logint) \
  shared(alpha, ast, psi) \
  schedule(static)
  for (i = 0; i < nbins; i++) {
    psi_i = psi[i];
    logc = (log(sin(alpha * psi_i)) - log(sin(psi_i))) / (1 - alpha) +
      log(sin((1 - alpha) * psi_i)) - log(sin(alpha * psi_i));
    logint = logc - exp(logc) * pow(ast, (- alpha / (1 - alpha)));
    integral[i] = exp(logint) * bin_width[i];
  }
  ll = llst + log(sum(integral));

  return ll;
}
