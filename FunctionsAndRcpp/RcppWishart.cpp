// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <iostream>
#include <vector>
using namespace Rcpp;
using namespace arma;

void showMatrix(mat X) {
  Rcout << "Armadillo matrix is" << std::endl << X << std::endl;
}

void showVector(rowvec x) { 
  Rcout << "Armadillo vector is " << std::endl << x << std::endl;
}


// [[Rcpp::export]]
extern "C" SEXP TransformObsWishart(SEXP observation, SEXP ind, SEXP mu0, SEXP sigma0){
  /* Function which performs the transformation of a new observation to eta_i,t in the thesis, for a given index.
   * 
   * Args:
   *     observation: new observation of the same length and order as mu0
   *     ind: which index is the transformation going to be performed for
   *     mu0: the in control mean
   *     sigma0: the in control covariance matrix. It is assumed that this is positive definite.
   *     
   * Returns:
   *     A vector which is the quantity eta_{index, t} in the thesis. 
   */

  typedef std::vector<double> stdvec;
  // inits
  int Ind = as<int>(ind)-1;
  colvec X_t = as<colvec>(observation);
  
  mat Sigma0 = as<mat>(sigma0);
  colvec Mu0 = as<colvec>(mu0);
  
  X_t = X_t - Mu0;
  // create matrices and create partitions of these.
  mat V_t = X_t*X_t.t();
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
  // Construct Schur complement
  mat Sigma_star = Sigma0-Sigma_i*Sigma_i.t()/sigma_ii;
  
  colvec eta_i = chol(Sigma_star.i()).t()*(V_ti/v_ii -Sigma_i/sigma_ii)*pow(v_ii,0.5);
  return wrap(eta_i);
}


  