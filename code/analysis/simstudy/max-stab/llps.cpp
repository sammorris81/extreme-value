// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

double integral_component(double *psi_i, double *alpha, double *ast_star) {
  double logc; double logint;
  double lsapsi_i = log(sin((*alpha) * (*psi_i)));

  logc = (lsapsi_i - log(sin(*psi_i))) / (1 - (*alpha)) +
          log(sin((1 - (*alpha)) * (*psi_i))) - lsapsi_i;

  logint = logc - exp(logc) * *ast_star;

  return logint;
}

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
  double integral; double llst; double logint;
  double ast; double psi_i;
  arma::mat ll(ns, nt); // arma::vec psi(nbins);
  double ast_star; 

  // psi = dscal(nbins, PI, mid_points, 1);
  // psi = PI * mid_points;

#pragma omp parallel for \
  private(t, i, ast, ast_star, llst, logint, integral, psi_i) \
  shared(a, alpha, ll, psi)
  for (s = 0; s < ns; s++) {
    for (t = 0; t < nt; t++) {
      ast = a(s, t);
      ast_star = pow(ast, (-alpha / (1 - alpha)));
      llst = log(alpha) - log(1 - alpha) - log(ast) / (1 - alpha);
      integral = 0;
      for (i = 0; i < nbins; i++) {
        psi_i = psi[i];
        logint = integral_component(&psi_i, &alpha, &ast_star);
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
  // arma::vec integral(nbins);
  double integral;
  double llst; double psi_i; double logint;
  double ll; double ast_star; 

  // psi = dscal(nbins, PI, mid_points, 1);
  ast_star = pow(ast, (-alpha / (1 - alpha)));
  llst = log(alpha) - log(1 - alpha) - log(ast) / (1 - alpha);

  integral = 0;
#pragma omp parallel for reduction(+:integral) \
  private(psi_i, logint) \
  shared(alpha, psi, ast_star)
  for (i = 0; i < nbins; i++) {
    psi_i = psi[i];
    logint = integral_component(&psi_i, &alpha, &ast_star);
    integral += exp(logint) * bin_width[i];
  }
  ll = llst + log(integral);

  return ll;
}
