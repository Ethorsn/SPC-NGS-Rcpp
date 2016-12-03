// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <iostream>
#include <math.h>       /* exp */

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
extern "C" SEXP cpd_likelihood_estimaton(SEXP data, SEXP mu0, SEXP Sigma0){
  /* Function which estimates the change point using a generalized likelihood ratio.
  * TODO: add bootstrapped confidence interval.
  * 
  * Args: 
  *  data: matrix containing the data from Phase 2 monitoring.
  *  mu0: in control mean.
  *  Sigma0: in control covariance matrix.
  *  
  * Returns:
  *  a vector containing values of the generalized likelihood ratio.
  */  
  
  rowvec mu = as<rowvec>(mu0);
  mat sigma = as<mat>(Sigma0);
  mat X = as<mat>(data);
  rowvec mu1(sigma.n_cols), ret(X.n_rows);
  double lik_val, c;
  for(int r=0; r<X.n_rows; r++){
    mat d = X.rows(r, X.n_rows-1);
    mu1 = mean(d);
    mat d1 = d.each_row()-mu1;
    mat d2 = d.each_row()-mu;
    c = -1*trace(d1*sigma.i()*d1.t() -
      d2*sigma.i()*d2.t())/2;
    lik_val = exp (c);
    ret.at(r) = lik_val;
  }
  return wrap(ret);
}

/***R
# library(mvnfast)
cpd_likelihood_estimaton(rbind(rmvn(20, mu=rep(0,3), sigma = diag(3)),rmvn(20, mu=rep(1,3), sigma = diag(3))),
                         mu0=rep(0,3),
                         Sigma0=diag(3))
*/

rowvec cpd_likelihood_estimaton_cpp(mat X, rowvec mu, mat sigma){
  /* Function which estimates the change point using a generalized likelihood ratio.
  * TODO: add bootstrapped confidence interval.
  * 
  * Args: 
  *  data: matrix containing the data from Phase 2 monitoring.
  *  mu0: in control mean.
  *  Sigma0: in control covariance matrix.
  *  
  * Returns:
  *  a vector containing values of the generalized likelihood ratio.
  */  
  rowvec mu1(sigma.n_cols), ret(X.n_rows);
  double lik_val, c;
  for(int r=0; r<X.n_rows; r++){
    mat d = X.rows(r, X.n_rows-1);
    mu1 = mean(d);
    mat d1 = d.each_row()-mu1;
    mat d2 = d.each_row()-mu;
    c = -1*trace(d1*sigma.i()*d1.t() -
      d2*sigma.i()*d2.t())/2;
    lik_val = exp (c);
    ret.at(r) = lik_val;
  }
  return ret;
}

mat mvrnormArma(int n, rowvec mu, mat sigma) {
  /* Function which simulates n multivariate normal obs using the cholesky decomposition.
  * 
  * Args:
  *     n: the number of observations to simulate.
  *     mu: mean of the normal distribution you want to simulate from.
  *     sigma: covariance matrix of the normal distribution you want to simulate from. It is assumed that it is positive definite.
  * 
  * Returns:
  *     n observations from a multivariate normal distribution with mean mu and covariance sigma.   
  */ 
  int ncols = sigma.n_cols;
  arma::mat Y = randn(n, ncols);
  mat Mu = repmat(mu, n,1);
  return Mu + Y * chol(sigma);
}

// [[Rcpp::export]]
extern "C" SEXP simulate_cpd_likelihood(SEXP NumbSims,SEXP N,SEXP tau, SEXP mu0, SEXP mu1, SEXP sigma0, SEXP sigma1, SEXP No_threads) {
  /* Function which simulates change points for different values of mu1 and sigma1.
  * 
  * Args:
  *     NumbSims: how many times to repeat the simulations
  *     N: how many OC observations to be simulated
  *     tau: how many observations before we go to OC 
  *     mu0: in control mean
  *     mu1: out-of-control mean 
  *     sigma0: in control covariance
  *     sigma1: out-of-control covariance
  *     
  *  Returns:
  *     A vector of simualted D, where D is the estimated change point minus the actual change point. 
  *     It represent how off our estimated changepoint is. 
  */
  int Nsims = as<int>(NumbSims);
  int n = as<int>(N);
  int TAU = as<int>(tau);
  rowvec Mu0 = as<rowvec>(mu0);
  rowvec Mu1 = as<rowvec>(mu1);
  mat Sigma0 = as<mat>(sigma0);
  mat Sigma1 = as<mat>(sigma1);
  int Nthread = as<int>(No_threads);
  
  rowvec ret_Vec = rowvec(Nsims);
  // specify matrices 
  int n_vars = Sigma0.n_cols;
  mat ICdata = mat(TAU,n_vars);
  mat OCdata = mat(n,n_vars);
  mat Alldata = mat(TAU+n,n_vars);
  size_t l;
//#pragma omp parallel for num_threads(Nthread) shared(ret_Vec) private(l,ICdata,OCdata,Alldata) firstprivate(TAU,n,Mu0,Mu1,Sigma0,Sigma1)
  for(l = 0; l < Nsims; l++){
    // simulate data
    ICdata = mvrnormArma(TAU, Mu0, Sigma0);
    OCdata = mvrnormArma(n, Mu1, Sigma1);
    // join and calculate changepoint
    Alldata = join_vert(ICdata, OCdata);
    rowvec ret_From_CDP = cpd_likelihood_estimaton_cpp(Alldata, Mu0, Sigma0);
    
    uword index;
    ret_From_CDP.max(index);
    int indexInt = index;
    
    ret_Vec[l] = indexInt-TAU;
  }
  return wrap(ret_Vec);
}

/***R
simulate_cpd_likelihood(NumbSims=10,N=20,tau=20,mu0=rep(0,4), mu1=rep(0.2,4), sigma0=diag(4), sigma1=diag(4), No_threads=2)
*/
