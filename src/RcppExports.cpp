// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Zdbcomparepairwise
NumericMatrix Zdbcomparepairwise(IntegerVector db, int nloci, bool display_progress = true);
RcppExport SEXP DNAprofiles_Zdbcomparepairwise(SEXP dbSEXP, SEXP nlociSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP );
        Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP );
        Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP );
        NumericMatrix __result = Zdbcomparepairwise(db, nloci, display_progress);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zdbcomparepairwisetrackhits
List Zdbcomparepairwisetrackhits(IntegerVector db, int nloci, int hit, bool display_progress = true);
RcppExport SEXP DNAprofiles_Zdbcomparepairwisetrackhits(SEXP dbSEXP, SEXP nlociSEXP, SEXP hitSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP );
        Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP );
        Rcpp::traits::input_parameter< int >::type hit(hitSEXP );
        Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP );
        List __result = Zdbcomparepairwisetrackhits(db, nloci, hit, display_progress);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zdbcomparepairwisemc
NumericMatrix Zdbcomparepairwisemc(IntegerVector db, int nloci, int njobs, int job);
RcppExport SEXP DNAprofiles_Zdbcomparepairwisemc(SEXP dbSEXP, SEXP nlociSEXP, SEXP njobsSEXP, SEXP jobSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP );
        Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP );
        Rcpp::traits::input_parameter< int >::type njobs(njobsSEXP );
        Rcpp::traits::input_parameter< int >::type job(jobSEXP );
        NumericMatrix __result = Zdbcomparepairwisemc(db, nloci, njobs, job);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zdbcomparepairwisemctrackhits
List Zdbcomparepairwisemctrackhits(IntegerVector db, int nloci, int hit, int njobs, int job);
RcppExport SEXP DNAprofiles_Zdbcomparepairwisemctrackhits(SEXP dbSEXP, SEXP nlociSEXP, SEXP hitSEXP, SEXP njobsSEXP, SEXP jobSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP );
        Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP );
        Rcpp::traits::input_parameter< int >::type hit(hitSEXP );
        Rcpp::traits::input_parameter< int >::type njobs(njobsSEXP );
        Rcpp::traits::input_parameter< int >::type job(jobSEXP );
        List __result = Zdbcomparepairwisemctrackhits(db, nloci, hit, njobs, job);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zdistapprox
List Zdistapprox(List dist, long maxn, double r0, double R, int method);
RcppExport SEXP DNAprofiles_Zdistapprox(SEXP distSEXP, SEXP maxnSEXP, SEXP r0SEXP, SEXP RSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type dist(distSEXP );
        Rcpp::traits::input_parameter< long >::type maxn(maxnSEXP );
        Rcpp::traits::input_parameter< double >::type r0(r0SEXP );
        Rcpp::traits::input_parameter< double >::type R(RSEXP );
        Rcpp::traits::input_parameter< int >::type method(methodSEXP );
        List __result = Zdistapprox(dist, maxn, r0, R, method);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zproductdist
List Zproductdist(NumericMatrix x, NumericMatrix prob, IntegerVector i, IntegerVector n, int N, double pr0, double prinf, bool returncumdist);
RcppExport SEXP DNAprofiles_Zproductdist(SEXP xSEXP, SEXP probSEXP, SEXP iSEXP, SEXP nSEXP, SEXP NSEXP, SEXP pr0SEXP, SEXP prinfSEXP, SEXP returncumdistSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type prob(probSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type n(nSEXP );
        Rcpp::traits::input_parameter< int >::type N(NSEXP );
        Rcpp::traits::input_parameter< double >::type pr0(pr0SEXP );
        Rcpp::traits::input_parameter< double >::type prinf(prinfSEXP );
        Rcpp::traits::input_parameter< bool >::type returncumdist(returncumdistSEXP );
        List __result = Zproductdist(x, prob, i, n, N, pr0, prinf, returncumdist);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ZfindIntervalcpp
IntegerVector ZfindIntervalcpp(NumericVector x, NumericVector breaks);
RcppExport SEXP DNAprofiles_ZfindIntervalcpp(SEXP xSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type breaks(breaksSEXP );
        IntegerVector __result = ZfindIntervalcpp(x, breaks);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ZcompKIpairswithtable
NumericVector ZcompKIpairswithtable(List X, IntegerMatrix db1, IntegerMatrix db2);
RcppExport SEXP DNAprofiles_ZcompKIpairswithtable(SEXP XSEXP, SEXP db1SEXP, SEXP db2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type X(XSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type db1(db1SEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type db2(db2SEXP );
        NumericVector __result = ZcompKIpairswithtable(X, db1, db2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ZcompKItargetsdbwithtable
NumericMatrix ZcompKItargetsdbwithtable(List X, IntegerMatrix db1, IntegerMatrix db2);
RcppExport SEXP DNAprofiles_ZcompKItargetsdbwithtable(SEXP XSEXP, SEXP db1SEXP, SEXP db2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type X(XSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type db1(db1SEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type db2(db2SEXP );
        NumericMatrix __result = ZcompKItargetsdbwithtable(X, db1, db2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ZcompKIwithtable
NumericVector ZcompKIwithtable(List X, IntegerMatrix db);
RcppExport SEXP DNAprofiles_ZcompKIwithtable(SEXP XSEXP, SEXP dbSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type X(XSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type db(dbSEXP );
        NumericVector __result = ZcompKIwithtable(X, db);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zrmpcpp
NumericVector Zrmpcpp(IntegerMatrix db, NumericMatrix f);
RcppExport SEXP DNAprofiles_Zrmpcpp(SEXP dbSEXP, SEXP fSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type db(dbSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type f(fSEXP );
        NumericVector __result = Zrmpcpp(db, f);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zstl_nth_element
NumericVector Zstl_nth_element(NumericVector x, int n);
RcppExport SEXP DNAprofiles_Zstl_nth_element(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        NumericVector __result = Zstl_nth_element(x, n);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Zsumprodxy
double Zsumprodxy(NumericVector x, NumericVector y);
RcppExport SEXP DNAprofiles_Zsumprodxy(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        double __result = Zsumprodxy(x, y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
