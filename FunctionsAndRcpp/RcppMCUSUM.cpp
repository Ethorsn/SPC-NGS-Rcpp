// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
extern "C" SEXP MCUSUM(SEXP Obs,SEXP Sold, SEXP mu0, SEXP Sigma0, SEXP k) {
  /* Function which computes Croisers MCUSUM scheme presented in the thesis.
   * 
   * Args:
   *     Obs: a new observation of the same length and order as mu0
   *     Sold: the vector S_{t} from previous computations
   *     mu0: the in control mean
   *     Sigma0: the in control covariance matrix
   *     k: the allowance constant
   *     
   * Returns:
   *     A list containing the charting statistic (named H_t in the thesis) and the new vector S_{t+1}
   */
  double K= as<double>(k);
  rowvec x = as<rowvec>(Sold);
  rowvec y = as<rowvec>(Obs);
  rowvec Mu0 = as<rowvec>(mu0);
  mat Sigma = as<mat>(Sigma0);
  double C;
  rowvec SNEW = zeros<rowvec>(Sigma.n_cols);
  
  mat SigmaInv = Sigma.i();
  y = y - Mu0;
  colvec Y = solve(Sigma , trans(x+y));
  // statistical distance, square root of Mahalanobis distance.
  double D = pow(as_scalar((x+y)*Y),0.5);
  
  if (D>K){
    double Dinv = as_scalar(pow(D,-1));
    SNEW = (x+y)*(1-K*Dinv);
    C = as_scalar(SNEW * SigmaInv * trans(SNEW));
  }else{
    C = 0;
  };
  return List::create(Named("Snew")=SNEW,
                      Named("Cstatistic")=C);
}

