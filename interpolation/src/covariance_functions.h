#ifndef COVARIANCE_FUNCTIONS_H
#define COVARIANCE_FUNCTIONS_H

#include <armadillo>
using namespace arma;

mat cov_Mat_cpp(const mat& ds, double sig2, double nugget,
                double phi, double kappa, bool diagonal = true);

mat cov_Sph_cpp(const mat& ds, double sig2, double nugget,
                double phi, bool diagonal = true);

mat cov_Exp_cpp(const mat& ds, double sig2, double nugget,
                double phi, bool diagonal = true);

mat cov_Gau_cpp(const mat& ds, double sig2, double nugget,
                double phi, bool diagonal = true);

#endif
