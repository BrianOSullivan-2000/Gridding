#include <RcppArmadillo.h>
#include "covariance_functions.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
mat cov_Mat_cpp(const mat& ds, double sig2, double nugget, 
                double phi, double kappa, bool diagonal) {
  
  int n = ds.n_rows;
  int m = ds.n_cols;
  mat V(n, m);
  double c = sig2 / (std::pow(2.0, kappa-1.0) * std::tgamma(kappa));
  
  if(diagonal && n == m) {
    for(int i=0; i<n; i++){
      V(i,i) = sig2 + nugget;
      for(int j=i+1; j<n; j++){
        double h = ds(i,j) / phi;
        double val = c * std::pow(h, kappa) * R::bessel_k(h, kappa, 1);
        V(i,j) = val;
        V(j,i) = val;
      }
    }
  } else {
    for(int i=0; i<n; i++){
      for(int j=0; j<m; j++){
        double h = ds(i,j) / phi;
        V(i,j) = c * std::pow(h, kappa) * R::bessel_k(h, kappa, 1);
      }
    }
  }
  return V;
}


// [[Rcpp::export]]
mat cov_Sph_cpp(const mat& ds, double sig2, double nugget, 
                double phi, bool diagonal) {
  
  int n = ds.n_rows;
  int m = ds.n_cols;
  mat V(n, m, fill::zeros);
  
  if (diagonal && n == m) {
    for (int i = 0; i < n; i++) {
      V(i,i) = sig2 + nugget;
      for (int j = i+1; j < n; j++) {
        double h = ds(i,j) / phi;
        if (h <= 1.0) {
          double h3 = h*h*h;
          double val = sig2 * (1.0 - 1.5*h + 0.5*h3);
          V(i,j) = val;
          V(j,i) = val;
        }
      }
    }
  } else {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        double h = ds(i,j) / phi;
        if (h <= 1.0) {
          double h3 = h*h*h;
          V(i,j) = sig2 * (1.0 - 1.5*h + 0.5*h3);
        }
      }
    }
  }
  return V;
}


// [[Rcpp::export]]
mat cov_Exp_cpp(const mat& ds, double sig2, double nugget, 
                double phi, bool diagonal) {
  
  int n = ds.n_rows;
  int m = ds.n_cols;
  mat V(n, m, fill::zeros);
  
  if (diagonal && n == m) {
    for (int i = 0; i < n; i++) {
      V(i,i) = sig2 + nugget;
      for (int j = i+1; j < n; j++) {
        double h = ds(i,j) / phi;
        double val = sig2 * std::exp(-h);
        V(i,j) = val;
        V(j,i) = val;
      }
    }
  } else {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        double h = ds(i,j) / phi;
        V(i,j) = sig2 * std::exp(-h);
      }
    }
  }
  return V;
}


// [[Rcpp::export]]
mat cov_Gau_cpp(const mat& ds, double sig2, double nugget, 
                double phi, bool diagonal) {
  
  int n = ds.n_rows;
  int m = ds.n_cols;
  mat V(n, m, fill::zeros);
  
  if (diagonal && n == m) {
    for (int i = 0; i < n; i++) {
      V(i,i) = sig2 + nugget;
      for (int j = i+1; j < n; j++) {
        double h = ds(i, j) / phi;
        double val = sig2 * std::exp(-h * h);
        V(i,j) = val;
        V(j,i) = val;
      }
    }
  } else {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        double h = ds(i, j) / phi;
        V(i, j) = sig2 * std::exp(-h * h);
      }
    }
  }
  
  // If nugget is exactly zero, add a small jitter so matrix is PSD
  if (nugget == 0.0) {
    double eps = 1e-12;
    V.diag() += eps;
  }
  
  return V;
}