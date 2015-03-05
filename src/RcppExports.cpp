// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// chisq_statistic
double chisq_statistic(std::vector<double> observed_left, std::vector<double> observed_right, size_t k, size_t n);
RcppExport SEXP maxchisq_chisq_statistic(SEXP observed_leftSEXP, SEXP observed_rightSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< std::vector<double> >::type observed_left(observed_leftSEXP );
        Rcpp::traits::input_parameter< std::vector<double> >::type observed_right(observed_rightSEXP );
        Rcpp::traits::input_parameter< size_t >::type k(kSEXP );
        Rcpp::traits::input_parameter< size_t >::type n(nSEXP );
        double __result = chisq_statistic(observed_left, observed_right, k, n);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pmaxchisq_permutation_internal
std::vector<double> pmaxchisq_permutation_internal(std::vector<double> b, size_t k, size_t n, size_t num_permutations, std::vector<size_t> y, std::vector<size_t> num_left, std::vector<size_t> class_counts);
RcppExport SEXP maxchisq_pmaxchisq_permutation_internal(SEXP bSEXP, SEXP kSEXP, SEXP nSEXP, SEXP num_permutationsSEXP, SEXP ySEXP, SEXP num_leftSEXP, SEXP class_countsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< std::vector<double> >::type b(bSEXP );
        Rcpp::traits::input_parameter< size_t >::type k(kSEXP );
        Rcpp::traits::input_parameter< size_t >::type n(nSEXP );
        Rcpp::traits::input_parameter< size_t >::type num_permutations(num_permutationsSEXP );
        Rcpp::traits::input_parameter< std::vector<size_t> >::type y(ySEXP );
        Rcpp::traits::input_parameter< std::vector<size_t> >::type num_left(num_leftSEXP );
        Rcpp::traits::input_parameter< std::vector<size_t> >::type class_counts(class_countsSEXP );
        std::vector<double> __result = pmaxchisq_permutation_internal(b, k, n, num_permutations, y, num_left, class_counts);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
