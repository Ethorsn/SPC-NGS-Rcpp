// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

void showMatrix(mat X) {
  Rcout << "Armadillo matrix is" << std::endl << X << std::endl;
}

void showVector(rowvec x) { 
  Rcout << "Armadillo vector is " << std::endl << x << std::endl;
}

// [[Rcpp::export]]
extern "C" SEXP CPDestimation(SEXP data, SEXP mu0, SEXP Sigma0){
  /* Function which estimates the change point using Hotellings two sample statistic.
   * 
   * Args:
   *     data: the data from Phase 2 monitoring which is to be used when estimating the change point.
   *     mu0: the in control mean
   *     Sigma0: the in control covariance
   *     
   * Returns:
   *     A vector where each value is Hotellings two sample statistic for each given change point. The position of 
   *     the largest value is the estimated change point.
   */
   
   
  mat Data = as<mat>(data);
  rowvec Mu0 = as<rowvec>(mu0);
  mat sigma0 = as<mat>(Sigma0);
  /* 
   * Flip the rows of the matrix Data. 
   * The most recent observation is then the first row etc.
   */
  Data = flipud(Data);
  mat Mean_mat = cumsum(Data,0);
  rowvec ret(Mean_mat.n_rows-1);
  /*
   * T = 1,2,3,...,N-1
   * if we flip data we T starts at one (second observation), 
   */
  int N = as_scalar(Mean_mat.n_rows);
  for(int T = 1; T <N; T++){
    rowvec y_T = (Mu0-Mean_mat.row(T)/(T+1))*(T*(N-T))/N;
    
    colvec Y = solve(sigma0,y_T.t());
    ret[T-1] = as_scalar(y_T*Y);
  };
  // flip back to the initial setting
  ret = fliplr(ret);
  return wrap(ret);
}


rowvec CPDestimation_cpp(mat data, rowvec mu0, mat Sigma0){
  /* Function which estimates the change point using Hotellings two sample statistic. USED IN C++ simulations.
   * 
   * Args:
   *     data: the data from Phase 2 monitoring which is to be used when estimating the change point.
   *     mu0: the in control mean
   *     Sigma0: the in control covariance
   *     
   * Returns:
   *     A vector where each value is Hotellings two sample statistic for each given change point. The position of 
   *     the largest value is the estimated change point.
   */
  
  /* 
   * Flip the rows of the matrix Data. 
   * The most recent observation is then the first row etc.
   */
  data = flipud(data);
  mat Mean_mat = cumsum(data,0);
  rowvec ret(Mean_mat.n_rows-1);
  /*
   * T = 1,2,3,...,N-1
   * if we flip data we T starts at one (second observation), 
   */
  int N = as_scalar(Mean_mat.n_rows);
  for(int T = 1; T <N; T++){
    rowvec y_T = (mu0-Mean_mat.row(T)/(T+1))*(T*(N-T))/N;
    
    colvec Y = solve(Sigma0,y_T.t());
    ret[T-1] = as_scalar(y_T*Y);
  };
  // flip back to the initial setting
  ret = fliplr(ret);
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
extern "C" SEXP SimulateCPDs(SEXP NumbSims,SEXP N,SEXP tau, SEXP mu0, SEXP mu1, SEXP sigma0, SEXP sigma1, SEXP No_threads) {
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
  int l;
#pragma omp parallel for num_threads(Nthread) shared(ret_Vec) private(l,ICdata,OCdata,Alldata) firstprivate(TAU,n,Mu0,Mu1,Sigma0,Sigma1)
  for(l = 0; l < Nsims; ++l){
    // simulate data
    ICdata = mvrnormArma(TAU, Mu0, Sigma0);
    OCdata = mvrnormArma(n, Mu1, Sigma1);
    // join and calculate changepoint
    Alldata = join_vert(ICdata, OCdata);
    rowvec ret_From_CDP = CPDestimation_cpp(Alldata, Mu0, Sigma0);
    
    uword index;
    ret_From_CDP.max(index);
    int indexInt = index;
    
    ret_Vec[l] = indexInt-TAU;
  }
  return wrap(ret_Vec);
}
/***R
#CPDestimation(matrix(c(1,2,3,1,2,3,1,2,3),byrow = FALSE, ncol=3),mu0=rep(0,3), Sigma0=diag(3)) # ok!
#SimulateCPDs(NumbSims = 10, N=70, tau=20, mu0=rep(0,3), mu1 = rep(10,3), sigma0 = diag(3), sigma1 = diag(3), No_threads = 2)
#SimulateCPDs(NumbSims = 1000, N=20, tau=20, mu0=rep(0,3), mu1 = rep(0.6,3), sigma0 = diag(3), sigma1 = diag(3))
*/
  