#include <RcppArmadillo.h>
#include "covariance_functions.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double spatial_likelihood_cpp(vec pars, mat ds, vec y, mat X, 
                              std::string cov_function, vec beta = vec()){
  
  // Parameters
  double sig2   = pars(0);
  double nugget = pars(1);
  double phi    = pars(2);
  
  // Get covariance matrix
  int n = ds.n_rows;
  mat V(n, n);
  
  if (cov_function == "Mat") {
    double kappa  = pars(3);
    V = cov_Mat_cpp(ds, sig2, nugget, phi, kappa, true);
    
  } else if (cov_function == "Sph") {
    V = cov_Sph_cpp(ds, sig2, nugget, phi, true);
    
  } else if (cov_function == "Exp") {
    V = cov_Exp_cpp(ds, sig2, nugget, phi, true);
    
  } else if (cov_function == "Gau") {
    V = cov_Gau_cpp(ds, sig2, nugget, phi, true);
    
  } else {
    Rcpp::stop("Unknown cov_function: %s", cov_function);
  }
  
    // Calculate likelihood
  mat U = chol(V);
  double logdetV = 2.0 * sum(log(U.diag()));
  
  // Get beta if not provided (GLS)
  if(beta.n_elem == 0){
    vec Vinv_y = solve(trimatu(U), solve(trimatl(U.t()), y));
    mat Vinv_X = solve(trimatu(U), solve(trimatl(U.t()), X));
    mat Xt_Vinv_X = X.t() * Vinv_X;
    vec Xt_Vinv_y = X.t() * Vinv_y;
    beta = solve(Xt_Vinv_X, Xt_Vinv_y);
  }
  
  vec res = y - X * beta;
  vec Vinv_res = solve(trimatu(U), solve(trimatl(U.t()), res));
  
  // Final calculation
  double quadform = dot(res, Vinv_res);
  double negloglik = 0.5 * (logdetV + quadform);

  return negloglik;
} 
