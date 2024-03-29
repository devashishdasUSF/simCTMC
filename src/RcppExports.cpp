// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// intensityFunc
double intensityFunc(double t, double baserate);
RcppExport SEXP _simCTMC_intensityFunc(SEXP tSEXP, SEXP baserateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type baserate(baserateSEXP);
    rcpp_result_gen = Rcpp::wrap(intensityFunc(t, baserate));
    return rcpp_result_gen;
END_RCPP
}
// Qmat
arma::mat Qmat(double t, double a, double b, double c);
RcppExport SEXP _simCTMC_Qmat(SEXP tSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(Qmat(t, a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// Sim_first_State
double Sim_first_State();
RcppExport SEXP _simCTMC_Sim_first_State() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(Sim_first_State());
    return rcpp_result_gen;
END_RCPP
}
// FiveStateSimulation
Rcpp::List FiveStateSimulation(double a, double b, double c);
RcppExport SEXP _simCTMC_FiveStateSimulation(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(FiveStateSimulation(a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// transition_times
Rcpp::List transition_times(double s0, double s1, double a, double b, double c);
RcppExport SEXP _simCTMC_transition_times(SEXP s0SEXP, SEXP s1SEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< double >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(transition_times(s0, s1, a, b, c));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_simCTMC_intensityFunc", (DL_FUNC) &_simCTMC_intensityFunc, 2},
    {"_simCTMC_Qmat", (DL_FUNC) &_simCTMC_Qmat, 4},
    {"_simCTMC_Sim_first_State", (DL_FUNC) &_simCTMC_Sim_first_State, 0},
    {"_simCTMC_FiveStateSimulation", (DL_FUNC) &_simCTMC_FiveStateSimulation, 3},
    {"_simCTMC_transition_times", (DL_FUNC) &_simCTMC_transition_times, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_simCTMC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
