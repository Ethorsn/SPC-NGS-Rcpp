// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;

rowvec SnewFun(rowvec Snew, rowvec Sold, rowvec mu0, mat Sigma0, double k) {
  /* Function which calculates the S_{t+1} from the MCUSUM scheme 
  * based on a Observation, S_old and the in control parameters and the allowance constant.  
  * 
  * Args:
  *     Observation: A vector of new observations of the same length and order as mu0
  *     Sold: A vector with values from the last iteration
  *     mu0: the in control mean vector
  *     Sigma0: the in control covariance matrix. Assumed to be positive definite
  *     k: the allowance constant 
  * 
  * Returns:
  *     A vector S which is the next iteration in the MCUSUM scheme.
  */
  //inits
  rowvec SNEW;
  
  colvec Y = solve(Sigma0, trans(Snew+Sold-mu0));
  // statistical distance, Mahalanobis distance.
  double D = pow(as_scalar((Snew+Sold-mu0)*Y),0.5);
  if (D>k){
    double Dinv = as_scalar(pow(D,-1));
    SNEW = (Snew+Sold-mu0)*(1-k*Dinv);
  }else{
    int length = Sigma0.n_rows;
    SNEW = zeros<rowvec>(length);
  };
  return SNEW;
}

double CFun(rowvec Snew, rowvec Sold, rowvec mu0, mat Sigma0, double k){
  /* Function which calculates the charting statistic from the MCUSUM scheme 
  * based on a Observation, S_old and the in control parameters and the allowance constant.  
  * 
  * Args:
  *     Observation: A vector of new observations of the same length and order as mu0
  *     Sold: A vector with values from the last iteration
  *     mu0: the in control mean vector
  *     Sigma0: the in control covariance matrix. Assumed to be positive definite
  *     k: the allowance constant 
  * 
  * Returns:
  *    A value of the charting statistic
  */
  //inits
  double C;
  rowvec SNEW;
  mat SigmaInv = Sigma0.i();
  
  colvec Y = solve(Sigma0, trans(Snew+Sold-mu0));
  // statistical distance, Mahalanobis distance.
  double D = pow(as_scalar((Snew+Sold-mu0)*Y),0.5);
  if (D>k){
    double Dinv = as_scalar(pow(D,-1));
    SNEW = (Snew+Sold-mu0)*(1-k*Dinv);
    C = as_scalar(SNEW * SigmaInv * trans(SNEW));
  }else{
    C = 0;
  };
  return C;
}
// [[Rcpp::export]]
extern "C" SEXP MCUSUM(SEXP Obs, SEXP Sold, SEXP mu0, SEXP Sigma0, SEXP k) {
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
  rowvec S_old = as<rowvec>(Sold);
  rowvec X = as<rowvec>(Obs);
  rowvec Mu0 = as<rowvec>(mu0);
  mat Sigma = as<mat>(Sigma0);
  double H;
  rowvec SNEW = zeros<rowvec>(Sigma.n_cols);
  
  mat SigmaInv = Sigma.i();
  X = X - Mu0;
  colvec Y = solve(Sigma , trans(X + S_old));
  // statistical distance, square root of Mahalanobis distance.
  double C = pow(as_scalar((X+S_old)*Y),0.5);
  if (C>K){
    double Cinv = as_scalar(pow(C,-1));
    SNEW = (X+S_old)*(1-K*Cinv);
    H = as_scalar(SNEW * SigmaInv * trans(SNEW));
  }else{
    H = 0;
  };
  return List::create(Named("Snew")=SNEW,
                      Named("Cstatistic")=H);
}