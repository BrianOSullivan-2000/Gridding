#include <RcppArmadillo.h>
#include "covariance_functions.h"
#include <chrono>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
vec predict_MVN_cpp(vec pars, mat ds, mat cross_ds,
                    vec y, mat X, mat new_X, 
                    std::string cov_function, vec beta = vec()){
  
  // Parameters
  double sig2   = pars(0);
  double nugget = pars(1);
  double phi    = pars(2);
  
  // Get covariance matrix
  int n = ds.n_rows;
  mat V(n, n);
  mat cross_V = cross_ds;
  
  if (cov_function == "Mat") {
    double kappa  = pars(3);
    V = cov_Mat_cpp(ds, sig2, nugget, phi, kappa, true);
    cross_V = cov_Mat_cpp(cross_ds, sig2, 0.0, phi, kappa, false);
    
  } else if (cov_function == "Sph") {
    V = cov_Sph_cpp(ds, sig2, nugget, phi, true);
    cross_V = cov_Sph_cpp(cross_ds, sig2, 0.0, phi, false);
    
  } else if (cov_function == "Exp") {
    V = cov_Exp_cpp(ds, sig2, nugget, phi, true);
    cross_V = cov_Exp_cpp(cross_ds, sig2, 0.0, phi, false);
    
  } else if (cov_function == "Gau") {
    V = cov_Gau_cpp(ds, sig2, nugget, phi, true);
    cross_V = cov_Gau_cpp(cross_ds, sig2, 0.0, phi, false);
    
  } else {
    Rcpp::stop("Unknown cov_function: %s", cov_function);
  }
  
  // Get inverse of covariance matrix
  mat U = chol(V);
  mat Vinv = solve(trimatu(U), solve(trimatl(U.t()), eye(n,n)));
  
  // Get beta if not provided (GLS)
  if(beta.n_elem == 0){
    vec Vinv_y = solve(trimatu(U), solve(trimatl(U.t()), y));
    mat Vinv_X = solve(trimatu(U), solve(trimatl(U.t()), X));
    mat Xt_Vinv_X = X.t() * Vinv_X;
    vec Xt_Vinv_y = X.t() * Vinv_y;
    beta = solve(Xt_Vinv_X, Xt_Vinv_y);
  }
  
  // Get residuals
  vec res = y - X * beta;
  
  // Predict at new locations
  vec pred = new_X * beta + cross_V.t() * (Vinv * res);
  
  return pred;
}
