// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
extern "C" SEXP HotellingT2(SEXP x_new, SEXP mu0, SEXP Sigma0) {
  /* Function which calculates Hotellings T2 statistic.
   * 
   * Args:
   *     x_new: A observation which should be on the same order as the data used to estimate mu0 and sigma0.
   *     mu0: The in control mean vector, estimated from a in control sample.
   *     Sigma0: The in control covariance matrix. It is assumed that this matrix is positive definite.
   * 
   * Returns:
   *     The value of Hotellings T2 statistic for this specific new observation.
   */
  
  // inits
  rowvec Xnew = as<rowvec>(x_new);
  rowvec Mu = as<rowvec>(mu0);
  mat Sigma = as<mat>(Sigma0);
  
  colvec Y = solve(Sigma, (Xnew-Mu).t());
  double ret = as_scalar((Xnew-Mu)*Y);
  return wrap(ret);
}


