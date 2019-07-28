#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>

// [[Rcpp::export]]
double intensityFunc(double t, double baserate = 1.0)
{
	return baserate*(2 - sin(2*M_PI*t));// +  sin(M_PI*(t-0.75)/2));
}