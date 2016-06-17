#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <vector>
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
    C = as_scalar(SNEW * SigmaInv * trans(SNEW));
  }else{
    C = 0;
  };
  return C;
}

rowvec TransformObsWishart(colvec Observation, int Ind, colvec Mu0, mat Sigma0){
  /* Function which performs the transformation of the new observation to eta_i,t in the thesis.
   * 
   * Args:
   *     Observation: A vector of new observations of the same length and order as mu0
   *     Ind: the index to perform the transformation on
   *     Mu0: in control mean
   *     Sigma0: in control covariance matrix
   *  
   * Returns:
   *     A vector which is the quantity eta_{i,t} for i=Ind. 
   */
  typedef std::vector<double> stdvec;
  
  Observation = Observation - Mu0;

  // create matrices and create partitions of these.
  mat V_t = Observation*Observation.t();
  // extract values at i'th diagonal
  double v_ii = V_t.at(Ind,Ind);
  double sigma_ii = Sigma0.at(Ind,Ind);
  
  // Extract the column and remove the element on the Ind'th position
  colvec Sigma_i = Sigma0.col(Ind);
  stdvec Sigma_i_tmp = conv_to< stdvec >::from(Sigma_i);
  // remove the ith element
  Sigma_i_tmp.erase(Sigma_i_tmp.begin()+Ind);
  /*
   *  We change the size of Sigma_i and dont really care which elements are left since 
   *  we are going to fill it with the results from Sigma_i_tmp.
   */
  Sigma_i.resize(Sigma_i.size()-1);
  Sigma_i = conv_to< colvec >::from(Sigma_i_tmp);
  
  // Do the same for the matrix V
  colvec V_ti = V_t.col(Ind);
  stdvec V_ti_tmp = conv_to< stdvec >::from(V_ti);
  V_ti_tmp.erase(V_ti_tmp.begin()+Ind);
  V_ti.resize(Sigma_i.size()-1);
  V_ti = conv_to< colvec >::from(V_ti_tmp);
  
  // remove Ind row and Ind col out of both matrices V and Sigma
  Sigma0.shed_col(Ind);
  Sigma0.shed_row(Ind);
  V_t.shed_col(Ind);
  V_t.shed_row(Ind);
  
  mat Sigma_star = Sigma0-Sigma_i*Sigma_i.t()/sigma_ii;
  
  colvec eta_i = chol(Sigma_star.t()).t()*(V_ti/v_ii - Sigma_i/sigma_ii)*pow(v_ii,0.5);
  return eta_i.t();
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
rowvec SimulateARL1Sigma(SEXP n, SEXP h, SEXP k, SEXP mu0, SEXP mu1, SEXP n0, SEXP Sigma0, SEXP Sigma1, SEXP No_threads) {
  /* Function which simulates the out-of-control control average run length (ARL1) or conditional expected delay for covariance under a given control limit, allowance constant
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
   *     Sigma1: the out-of-control covariance matrix. It is assumed that this is positive definite.
   *     No_threads: The number of threads to use in the simulations.
   * 
   * Returns:
   *     Returns a vector of length n which contains simulated values of the out-of-control average run length or conditional expected delay for the covariance.
   */

  // inits
  int N = as<int>(n);
  double H = as<double>(h);
  double K = as<double>(k);
  int Nthread = as<int>(No_threads);
  int N0 = as<int>(n0);
  rowvec Mu = as<rowvec>(mu0);
  rowvec Mu1 = as<rowvec>(mu1);
  
  mat Sigma = as<mat>(Sigma0);
  mat sigma1 = as<mat>(Sigma1);
  rowvec VecReturn = zeros<rowvec>(N);
  
  size_t l; 
  int counter, IamMax;
  double Cstat;
  rowvec S_old(Sigma.n_cols-1), Sims(Sigma.n_cols), TransformObs(Sigma.n_cols-1), TmpVector(Sigma.n_cols);
  
  mat UnitMatrix = eye(Sigma.n_cols-1,Sigma.n_cols-1);
  rowvec ZeroMean = zeros<rowvec>(Sigma.n_cols-1);
 
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
                                                       
#pragma omp parallel for num_threads(Nthread) shared(VecReturn) private(l, counter, Cstat, S_old, Sims,TransformObs, IamMax) firstprivate(H, K, Mu, Mu1, Sigma, sigma1, N0,UnitMatrix,ZeroMean)
      for (l = 0; l < N;++l){
        // Inits for a run
        counter = 1;
        Cstat = 0;
        S_old = zeros<rowvec>(Sigma.n_cols-1);
        // simulating the expectation.
        for (int k=0; k<10000; ++k){
          if (k<N0){
            Sims = mvrnormArma(Mu, Sigma);
          }
          else
          {
            Sims = mvrnormArma(Mu1, sigma1);
          }
          
          // Calculate all different Cstat for different indexes. Place these in TmpVector 
          for(int m=0; m < Sigma.n_cols; ++m){
            // Transform observation according to singular wishart properties
            TransformObs = TransformObsWishart(Sims.t(), m, Mu.t(), Sigma);
            
            // Cstat calculations
            Cstat = CFun(TransformObs, S_old, ZeroMean, UnitMatrix, K);
            // store Cstat
            TmpVector[m] = Cstat;
          };
          // Which value is the largest?
          uword index;
          Cstat = TmpVector.max(index);
          int IndexMe = index; 
          // Update S_old.
          TransformObs = TransformObsWishart(Sims.t(), IndexMe, Mu.t(), Sigma);
          S_old = SnewFun(TransformObs , S_old, ZeroMean, UnitMatrix, K);
          // did we break
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
    throw e;
    rowvec ret;
    return ret;
  }
  return VecReturn;
}
/***R
#SimulateARL1Sigma(n=2, h=10, k=0.1, mu0=rep(0,48), mu1=rep(0,48), n0=0, Sigma0=diag(48), Sigma1=diag(48), No_threads=3) 
#SimulateARL1Sigma(n=2, h=10, k=0.1, mu0=rep(0,48), mu1=rep(0,48), n0=0, Sigma0=diag(48), Sigma1=diag(48), No_threads=2) 
#(SEXP n, SEXP h, SEXP k, SEXP mu0, SEXP mu1, SEXP n0, SEXP Sigma0, SEXP Sigma1, SEXP No_threads)
*/