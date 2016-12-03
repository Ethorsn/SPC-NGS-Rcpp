#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

//#include "McusumHelper.h"
// error in dyn.load(), symbol not found, expected in: flat namespace, seem to be related to Rcpp?

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

rowvec SnewFun(rowvec Snew, rowvec Sold, rowvec mu0, mat Sigma0, double k) {
  /* Function which calculates the S_{t+1} from the MCUSUM control chart
  * based on a Snew, Sold, the in control parameters and the allowance constant.  
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
  mat SigmaInv = Sigma0.i();
  
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

rowvec mvrnormArma(rowvec mu, mat sigma) {
  /* Function which simulates ONE multivariate normal obs using the cholesky decomposition.
   * 
   * Args:
   *     mu: mean of the normal distribution you want to simulate from.
   *     sigma: covariance matrix of the normal distribution you want to simulate from. It is assumed that it is positive definite.
   * 
   * Returns:
   *     A observation from a multivariate normal distribution with mean mu and covariance sigma.   
   */ 
  
  int ncols = sigma.n_cols;
  rowvec Y = randn<rowvec>(ncols);
  rowvec ret = mu + Y * chol(sigma);
  return ret;
}

// Simulate ARL0 based on k and h in parallel using open MP. 
// [[Rcpp::export]]
rowvec SimulateARL1mean(SEXP n, SEXP h, SEXP k, SEXP mu0, SEXP mu1, SEXP n0, SEXP Sigma0, SEXP No_threads) {
  /* Function which simulates the out-of-control control average run length (ARL1) 
   * or conditional expected delay under a given control limit, allowance constant
   * and in control parameters.
   * 
   * Args:
   *     n: The number of simulations to perform.
   *     h: the control limit
   *     k: the allowance constant
   *     mu0: the in control mean
   *     mu1: the out-of-control mean
   *     n0: the number of observations before the process goes out-of-control
   *     Sigma0: the in control covariance matrix. It is assumed that this is positive definite.
   *     No_threads: The number of threads to use in the simulations.
   * 
   * Returns:
   *     Returns a vector of length n which contains simulated values of the 
   *     out-of-control average run length or conditional expected delay.
   */
  
  // inits
  int N = as<int>(n);
  double H = as<double>(h);
  double K = as<double>(k);
  int Nthread = as<int>(No_threads);
  int N0 = as<int>(n0);
  rowvec Mu = as<rowvec>(mu0);
  
  Mu.print("IC mean vector");
  rowvec Mu1 = as<rowvec>(mu1);
  Mu1.print("OC mean vector");
  
  mat Sigma = as<mat>(Sigma0);
  rowvec VecReturn = zeros<rowvec>(N);
  
  size_t l; 
  int counter;
  double Cstat;
  rowvec S_old(Sigma.n_cols), Sims(Sigma.n_cols);
  
  class InvalidDim_exception: public std::exception {};
  class NegParams_exception: public std::exception {};
                                                     
  try{
    if (Sigma.n_cols != Mu.size() || Sigma.n_rows != Mu.size() || Sigma.n_cols != Mu1.size() || Sigma.n_rows != Mu1.size())
    {
      // Sigma, Mu and Mu1 does not have the same dimensions.
      throw InvalidDim_exception();
    }
    else if (K < 0 || H < 0)
    {
      // k or h is negative
      throw NegParams_exception();
    }
    else
    {
      Rcout << "Inits done! Starting simulations. \n" << endl;
    }
    
    /*
     The addition of #pragma from OpenMP fixes parallel stuff...
     So the following code is in parallel, using Nthread (num_threads(Nthread)) threads.
     */
    
#pragma omp parallel for num_threads(Nthread) shared(VecReturn) private(l, counter, Cstat, S_old, Sims) firstprivate(H, K, Mu, Mu1, Sigma, N0)
    for (l = 0; l < N;++l){
      // Inits for a run
      counter = 1;
      Cstat = 0;
      S_old = zeros<rowvec>(Sigma.n_cols);
      // simulating the expectation.
      for (int k=0; k<10000; ++k){
        if (k<N0){
          Sims = mvrnormArma(Mu, Sigma);
        }
        else{
          Sims = mvrnormArma(Mu1, Sigma);
        }
        // Use MCUSUM scheme, first update Cstat then S_old.
        Cstat = CFun(Sims, S_old, Mu, Sigma, K);
        S_old = SnewFun(Sims, S_old, Mu, Sigma, K);
        if (Cstat > H){break;};
        counter += 1;
      };
      VecReturn(l) = counter-N0;    
    };
  }
  catch(NegParams_exception e){
    // add dimensions or print values of parameters
    cerr << "Negative values of the control limit (h) or allowance constant (k) are not allowed." << endl;
    throw e;
    rowvec ret;
    return ret;
  }
  catch(InvalidDim_exception e){
    // add dimensions or print values of parameters
    cerr << "The dimensions of the mean vector (mu0 or mu1) and Covariance matrix are not the same." << endl;
    rowvec ret;
    return ret;
  }
  return VecReturn;
}

