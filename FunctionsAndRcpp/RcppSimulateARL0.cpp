#include <RcppArmadillo.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
rowvec SnewFun(rowvec Observation, rowvec Sold, rowvec mu0, mat Sigma0, double k) {
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
  mat SigmaInv = Sigma0.i();
  
  colvec Y = solve(Sigma0, trans(Observation+Sold-mu0));
  // statistical distance, Mahalanobis distance.
  double D = pow(as_scalar((Observation+Sold-mu0)*Y),0.5);
  if (D>k){
    double Dinv = as_scalar(pow(D,-1));
    SNEW = (Observation+Sold-mu0)*(1-k*Dinv);
  }else{
    int length = Sigma0.n_rows;
    SNEW = zeros<rowvec>(length);
  };
  return SNEW;
}

double CFun(rowvec Observation, rowvec Sold, rowvec mu0, mat Sigma0, double k){
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
  
  colvec Y = solve(Sigma0, trans(Observation+Sold-mu0));
  // statistical distance, Mahalanobis distance.
  double D = pow(as_scalar((Observation+Sold-mu0)*Y),0.5);
  if (D>k){
    double Dinv = as_scalar(pow(D,-1));
    SNEW = (Observation+Sold-mu0)*(1-k*Dinv);
    C = as_scalar(SNEW * SigmaInv * SNEW.t());
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
rowvec SimulateARL0(SEXP n, SEXP h, SEXP k, SEXP mu0, SEXP Sigma0, SEXP No_threads) {
  /* Function which simulates the in control average run length (ARL0) under a given control limit, allowance constant
   * and in control parameters.
   * 
   * Args:
   *     n: The number of simulations to perform.
   *     h: the control limit
   *     k: the allowance constant
   *     mu0: the in control mean
   *     Sigma0: the in control covariance matrix. It is assumed that this is positive definite.
   *     No_threads: The number of threads to use in the simulations.
   * 
   * Returns:
   *     Returns a vector of length n which contains simulated values of the average run length.
   */
  
  // inits
  int N = as<int>(n);
  double H = as<double>(h);
  double K = as<double>(k);
  int Nthread = as<int>(No_threads);
  
  rowvec Mu = as<rowvec>(mu0);
  
  mat Sigma = as<mat>(Sigma0);
  rowvec VecReturn = zeros<rowvec>(N);
  
  size_t l; 
  int counter;
  double Cstat;
  rowvec S_old(Sigma.n_cols), Sims(Sigma.n_cols);
  // Define exceptions.
  class InvalidDim_exception: public std::exception {};
  class NegParams_exception: public std::exception {};
                                                     
  try{
    if (Sigma.n_cols != Mu.size() || Sigma.n_rows != Mu.size())
    {
      // Sigma, Mu and does not have the same dimensions.
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
     So the following code is in parallel, using num_threads(Nthread) threads.
     */
#if _OPENMP
#pragma omp parallel for num_threads(Nthread) shared(VecReturn) private(l, counter, Cstat, S_old, Sims) firstprivate(H, K, Mu, Sigma)
#endif   
    for (l = 0; l < N;++l){
      // Inits for a run
      counter = 1;
      Cstat = 0;
      S_old = zeros<rowvec>(Sigma.n_cols);
      // simulating the expectation.
      for (int k=0; k<10000; ++k){
        Sims = mvrnormArma(Mu, Sigma);
        // Use MCUSUM scheme, first update Cstat then S_old.
        Cstat = CFun(Sims, S_old, Mu, Sigma, K);
        S_old = SnewFun(Sims, S_old, Mu, Sigma, K);
        if (Cstat > H){break;};
        counter += 1;
      };
      VecReturn(l) = counter;  
    };
  }
  catch(NegParams_exception e){
    // Negative values
    Rcout << "Negative values of the control limit (h) or allowance constant (k) are not allowed." << endl;
    throw e;
    rowvec ret;
    return ret;
  }
  catch(InvalidDim_exception e){
    // Dimensions or print values of parameters
    Rcout << "The dimensions of the mean vector (mu0 or mu1) and Covariance matrix are not the same." << endl;
    throw e;
    rowvec ret;
    return ret;
  }
  return VecReturn;
}

/*** R
#SimulateARL0(n=10,h=5, k=1, mu0 = rep(0,5), Sigma0=diag(2),1) # error 
#SimulateARL0(n=10,h=-1, k=1, mu0 = rep(0,5), Sigma0=diag(5),1) # error
#SimulateARL0(n=10,h=20, k=0.3, mu0 = rep(0,5), Sigma0=diag(5),3) # error
# Works!
#SimulateARL0(n=100,h=0.01, k=0.3, mu0 = rep(0,5), Sigma0=diag(5),3)
*/
