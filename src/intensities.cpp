#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>

// [[Rcpp::export]]
double intensityFunc(double t, double baserate = 1.0)
{
	return baserate*(2 - sin(2*M_PI*t));// +  sin(M_PI*(t-0.75)/2));
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <iterator>
#include <algorithm>


using namespace Rcpp;


double q_12(double x, double a = 1., double b = 0., double c = 2.) {return a*(sin(-2*M_PI*(x-b)) + c);}
double q_13(double x) {return cos(2.0*M_PI*x) + 1.0;}
double q_21 (double x) {return sin(3.0*M_PI*x) + 1.0;}
double q_23 (double x) {return sin(2.0*M_PI*(x-0.1/5)) + 1.0;}
double q_31 (double x) {return 2.*sin(2.0*M_PI*x) + 2.;}
double q_32 (double x) {return 2.*cos(2.*M_PI*(x+.1/4)) + 2.;}

// [[Rcpp::export]] 
arma::mat Qmat(double t, double a = 1., double b = 0., double c = 2.0) {
  arma::mat Q(5,5);
  Q.zeros();
  Q(0,0) = -(q_12(t,a,b,c)+q_13(t));
  Q(0,1) = q_12(t,a,b,c); Q(0,2) = q_13(t);
  Q(1,0) = q_21(t); Q(1,1) = -(q_21(t)+q_23(t)); Q(1,2) = q_23(t);
  Q(2,0) = q_31(t); Q(2,1) = q_32(t);
  Q(2,2) = -2.*(q_31(t)+q_32(t)); Q(2,3) = q_32(t); Q(2,4) = q_31(t);
  Q(3,2) = q_12(t,a,b,c); Q(3,3) = -(q_12(t,a,b,c)+q_13(t)); Q(3,4) = q_13(t);
  Q(4,2) = q_21(t); Q(4,3) = q_23(t); Q(4,4) = -(q_21(t)+q_23(t));
  return Q;
}

const double lambda = 100.0;
const Rcpp::NumericVector StateSpace = Rcpp::NumericVector::create(1., 2., 3., 4., 5.);
// const Rcpp::NumericVector initial_prob(5, 0.20000);

Rcpp::NumericVector Jump_times() {
  
  Rcpp::NumericVector A(1000);
  // double lambda = 100;
  A = Rcpp::rexp(1000, 100);
  
  Rcpp::NumericVector temp(1000);
  std::partial_sum(A.begin(), A.end(), temp.begin());
  
  return temp[temp <= 1.0];
}
// [[Rcpp::export]] 
double Sim_first_State() {
  Rcpp::NumericVector initial_prob(5, 0.20000);
  Rcpp::NumericVector First = Rcpp::RcppArmadillo::sample(StateSpace, 1, true, initial_prob);
  return First[0];
}

Rcpp::NumericVector Sim_next_State(double i, double JumpT, double a = 1., double b = 0., double c = 1.0) {
  arma::mat Pmat = arma::eye(5,5) + (1/lambda)*Qmat(JumpT,a,b,c);
  int j = floor(i);
  return Rcpp::RcppArmadillo::sample(StateSpace, 1, true, Rcpp::wrap(Pmat.row(j-1)));
}

// [[Rcpp::export]] 
Rcpp::List FiveStateSimulation(double a = 1., double b = 0., double c = 2.0) {
  Rcpp::NumericVector TimePoints = Jump_times();
  int N = TimePoints.size();
  double First = Sim_first_State();
  Rcpp::NumericVector States(N+1, First);
  
  if (N >= 1) {
    for(int i = 0; i< N; i ++) {
      States[i+1] = (Sim_next_State(States[i], TimePoints[i],a,b,c))[0];
    }
  }
  Rcpp::List result;
  result["States"] = States;
  result["JumpTimes"] = TimePoints;
  return result;
}

// [[Rcpp::export]]
Rcpp::List transition_times(double s0 = 1.0,
  double s1 = 2.0,
  double a = 1., 
  double b = 0., 
  double c = 2.0) {

  Rcpp::List A;
  A = FiveStateSimulation(a,b,c);
  Rcpp::NumericVector S = A[0];
  Rcpp::NumericVector T = A[1];
  Rcpp::NumericVector index(S.size()-1,0);
  Rcpp::NumericVector T0(S.size(),0);

  for(int i = 1; i< S.size(); i++) {
    T0[i] = T[i-1];
  }
  Rcpp::NumericVector Enter_s0(S.size(),0);
  Rcpp::NumericVector Exit_s0(S.size()-1,0);
  for(int i=0; i<S.size()-1; i++){
    if(S[i] == s0 && S[i+1] == s1) {
      index[i] = 1.0;
    }
    if(S[i] == s0 && S[i+1] != s0) {
      Exit_s0[i] = 1.0;
    }
    if(S[i] != s0 && S[i+1] == s0) {
      Enter_s0[i+1] = 1.0;
    }
  }
  if(S[0] == s0) Enter_s0[0] = 1.0;

  Rcpp::List result;

  result["time_points"] = T[index == 1];
  result["Enter_s0_points"] = T0[Enter_s0 == 1];
  result["Exit_s0_points"] = T[Exit_s0 == 1];

  return result;
}




