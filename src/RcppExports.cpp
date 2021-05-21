// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// anccall_c
void anccall_c(double deltaf, int window, int SMAX, std::string anclikdir, std::string output);
RcppExport SEXP _rLCLAE_anccall_c(SEXP deltafSEXP, SEXP windowSEXP, SEXP SMAXSEXP, SEXP anclikdirSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type deltaf(deltafSEXP);
    Rcpp::traits::input_parameter< int >::type window(windowSEXP);
    Rcpp::traits::input_parameter< int >::type SMAX(SMAXSEXP);
    Rcpp::traits::input_parameter< std::string >::type anclikdir(anclikdirSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    anccall_c(deltaf, window, SMAX, anclikdir, output);
    return R_NilValue;
END_RCPP
}
// ancfreq_c
void ancfreq_c(int n, int type, std::string pop1, std::string pop2, std::string input, std::string output);
RcppExport SEXP _rLCLAE_ancfreq_c(SEXP nSEXP, SEXP typeSEXP, SEXP pop1SEXP, SEXP pop2SEXP, SEXP inputSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< std::string >::type pop1(pop1SEXP);
    Rcpp::traits::input_parameter< std::string >::type pop2(pop2SEXP);
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    ancfreq_c(n, type, pop1, pop2, input, output);
    return R_NilValue;
END_RCPP
}
// ancfreq_merge
void ancfreq_merge(std::string hapfreq, std::string dipfreq, std::string outputdir, int type);
RcppExport SEXP _rLCLAE_ancfreq_merge(SEXP hapfreqSEXP, SEXP dipfreqSEXP, SEXP outputdirSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type hapfreq(hapfreqSEXP);
    Rcpp::traits::input_parameter< std::string >::type dipfreq(dipfreqSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputdir(outputdirSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    ancfreq_merge(hapfreq, dipfreq, outputdir, type);
    return R_NilValue;
END_RCPP
}
// anclik_c
void anclik_c(int sum, int type, int testid, std::string genolik, std::string ancfreq, std::string output);
RcppExport SEXP _rLCLAE_anclik_c(SEXP sumSEXP, SEXP typeSEXP, SEXP testidSEXP, SEXP genolikSEXP, SEXP ancfreqSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type sum(sumSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type testid(testidSEXP);
    Rcpp::traits::input_parameter< std::string >::type genolik(genolikSEXP);
    Rcpp::traits::input_parameter< std::string >::type ancfreq(ancfreqSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    anclik_c(sum, type, testid, genolik, ancfreq, output);
    return R_NilValue;
END_RCPP
}
// filt1_dip
void filt1_dip(int n, int skip, int dur, std::string input, std::string output);
RcppExport SEXP _rLCLAE_filt1_dip(SEXP nSEXP, SEXP skipSEXP, SEXP durSEXP, SEXP inputSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type dur(durSEXP);
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    filt1_dip(n, skip, dur, input, output);
    return R_NilValue;
END_RCPP
}
// filt1_hap
void filt1_hap(int n, int skip, int dur, std::string input, std::string output);
RcppExport SEXP _rLCLAE_filt1_hap(SEXP nSEXP, SEXP skipSEXP, SEXP durSEXP, SEXP inputSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< int >::type dur(durSEXP);
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    filt1_hap(n, skip, dur, input, output);
    return R_NilValue;
END_RCPP
}
// glpow
double glpow(double x);
RcppExport SEXP _rLCLAE_glpow(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(glpow(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rLCLAE_anccall_c", (DL_FUNC) &_rLCLAE_anccall_c, 5},
    {"_rLCLAE_ancfreq_c", (DL_FUNC) &_rLCLAE_ancfreq_c, 6},
    {"_rLCLAE_ancfreq_merge", (DL_FUNC) &_rLCLAE_ancfreq_merge, 4},
    {"_rLCLAE_anclik_c", (DL_FUNC) &_rLCLAE_anclik_c, 6},
    {"_rLCLAE_filt1_dip", (DL_FUNC) &_rLCLAE_filt1_dip, 5},
    {"_rLCLAE_filt1_hap", (DL_FUNC) &_rLCLAE_filt1_hap, 5},
    {"_rLCLAE_glpow", (DL_FUNC) &_rLCLAE_glpow, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rLCLAE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
