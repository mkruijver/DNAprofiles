// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Zprnextallele
NumericVector Zprnextallele(IntegerVector i, IntegerMatrix seen, NumericVector fr, double theta);
RcppExport SEXP _DNAprofiles_Zprnextallele(SEXP iSEXP, SEXP seenSEXP, SEXP frSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type seen(seenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fr(frSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Zprnextallele(i, seen, fr, theta));
    return rcpp_result_gen;
END_RCPP
}
// Zprnextalleles
NumericVector Zprnextalleles(IntegerMatrix ij, IntegerMatrix seen, NumericVector fr, double theta);
RcppExport SEXP _DNAprofiles_Zprnextalleles(SEXP ijSEXP, SEXP seenSEXP, SEXP frSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type ij(ijSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type seen(seenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fr(frSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Zprnextalleles(ij, seen, fr, theta));
    return rcpp_result_gen;
END_RCPP
}
// Zallfinitepos
bool Zallfinitepos(NumericVector x, int i1, int i2);
RcppExport SEXP _DNAprofiles_Zallfinitepos(SEXP xSEXP, SEXP i1SEXP, SEXP i2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type i1(i1SEXP);
    Rcpp::traits::input_parameter< int >::type i2(i2SEXP);
    rcpp_result_gen = Rcpp::wrap(Zallfinitepos(x, i1, i2));
    return rcpp_result_gen;
END_RCPP
}
// Zcountallelescpp
NumericMatrix Zcountallelescpp(IntegerMatrix x, int Amax, NumericVector w);
RcppExport SEXP _DNAprofiles_Zcountallelescpp(SEXP xSEXP, SEXP AmaxSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type Amax(AmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(Zcountallelescpp(x, Amax, w));
    return rcpp_result_gen;
END_RCPP
}
// Zdbcomparepairwise
NumericMatrix Zdbcomparepairwise(IntegerVector db, int nloci, bool display_progress);
RcppExport SEXP _DNAprofiles_Zdbcomparepairwise(SEXP dbSEXP, SEXP nlociSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP);
    Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(Zdbcomparepairwise(db, nloci, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// Zdbcomparepairwisetrackhits
List Zdbcomparepairwisetrackhits(IntegerVector db, int nloci, int hit, bool display_progress);
RcppExport SEXP _DNAprofiles_Zdbcomparepairwisetrackhits(SEXP dbSEXP, SEXP nlociSEXP, SEXP hitSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP);
    Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP);
    Rcpp::traits::input_parameter< int >::type hit(hitSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(Zdbcomparepairwisetrackhits(db, nloci, hit, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// Zdbcomparepairwisemc
NumericMatrix Zdbcomparepairwisemc(IntegerVector db, int nloci, int njobs, int job);
RcppExport SEXP _DNAprofiles_Zdbcomparepairwisemc(SEXP dbSEXP, SEXP nlociSEXP, SEXP njobsSEXP, SEXP jobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP);
    Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP);
    Rcpp::traits::input_parameter< int >::type njobs(njobsSEXP);
    Rcpp::traits::input_parameter< int >::type job(jobSEXP);
    rcpp_result_gen = Rcpp::wrap(Zdbcomparepairwisemc(db, nloci, njobs, job));
    return rcpp_result_gen;
END_RCPP
}
// Zdbcomparepairwisemctrackhits
List Zdbcomparepairwisemctrackhits(IntegerVector db, int nloci, int hit, int njobs, int job);
RcppExport SEXP _DNAprofiles_Zdbcomparepairwisemctrackhits(SEXP dbSEXP, SEXP nlociSEXP, SEXP hitSEXP, SEXP njobsSEXP, SEXP jobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type db(dbSEXP);
    Rcpp::traits::input_parameter< int >::type nloci(nlociSEXP);
    Rcpp::traits::input_parameter< int >::type hit(hitSEXP);
    Rcpp::traits::input_parameter< int >::type njobs(njobsSEXP);
    Rcpp::traits::input_parameter< int >::type job(jobSEXP);
    rcpp_result_gen = Rcpp::wrap(Zdbcomparepairwisemctrackhits(db, nloci, hit, njobs, job));
    return rcpp_result_gen;
END_RCPP
}
// Zdistapprox
List Zdistapprox(List dist, long maxn, double r0, double R, int method);
RcppExport SEXP _DNAprofiles_Zdistapprox(SEXP distSEXP, SEXP maxnSEXP, SEXP r0SEXP, SEXP RSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dist(distSEXP);
    Rcpp::traits::input_parameter< long >::type maxn(maxnSEXP);
    Rcpp::traits::input_parameter< double >::type r0(r0SEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(Zdistapprox(dist, maxn, r0, R, method));
    return rcpp_result_gen;
END_RCPP
}
// Zproductdist
List Zproductdist(NumericMatrix x, NumericMatrix prob, IntegerVector i, IntegerVector n, int N, double pr0, double prinf, bool returncumdist);
RcppExport SEXP _DNAprofiles_Zproductdist(SEXP xSEXP, SEXP probSEXP, SEXP iSEXP, SEXP nSEXP, SEXP NSEXP, SEXP pr0SEXP, SEXP prinfSEXP, SEXP returncumdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type prob(probSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type pr0(pr0SEXP);
    Rcpp::traits::input_parameter< double >::type prinf(prinfSEXP);
    Rcpp::traits::input_parameter< bool >::type returncumdist(returncumdistSEXP);
    rcpp_result_gen = Rcpp::wrap(Zproductdist(x, prob, i, n, N, pr0, prinf, returncumdist));
    return rcpp_result_gen;
END_RCPP
}
// Zexactq
double Zexactq(double t, NumericMatrix x, NumericMatrix prob, IntegerVector i, IntegerVector n, double pr0);
RcppExport SEXP _DNAprofiles_Zexactq(SEXP tSEXP, SEXP xSEXP, SEXP probSEXP, SEXP iSEXP, SEXP nSEXP, SEXP pr0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type prob(probSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type pr0(pr0SEXP);
    rcpp_result_gen = Rcpp::wrap(Zexactq(t, x, prob, i, n, pr0));
    return rcpp_result_gen;
END_RCPP
}
// ZfindIntervalcpp
IntegerVector ZfindIntervalcpp(NumericVector x, NumericVector breaks);
RcppExport SEXP _DNAprofiles_ZfindIntervalcpp(SEXP xSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(ZfindIntervalcpp(x, breaks));
    return rcpp_result_gen;
END_RCPP
}
// Zki
NumericMatrix Zki(IntegerMatrix x1, int manytomany, IntegerMatrix x2, IntegerVector x1ind, IntegerVector x2ind, NumericMatrix fr, double k0, double k1, double k2, double theta, int retpermarker);
RcppExport SEXP _DNAprofiles_Zki(SEXP x1SEXP, SEXP manytomanySEXP, SEXP x2SEXP, SEXP x1indSEXP, SEXP x2indSEXP, SEXP frSEXP, SEXP k0SEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP thetaSEXP, SEXP retpermarkerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< int >::type manytomany(manytomanySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x1ind(x1indSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x2ind(x2indSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type fr(frSEXP);
    Rcpp::traits::input_parameter< double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< double >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< double >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type retpermarker(retpermarkerSEXP);
    rcpp_result_gen = Rcpp::wrap(Zki(x1, manytomany, x2, x1ind, x2ind, fr, k0, k1, k2, theta, retpermarker));
    return rcpp_result_gen;
END_RCPP
}
// ZcompKIpairswithtable
NumericVector ZcompKIpairswithtable(List X, IntegerMatrix db1, IntegerMatrix db2);
RcppExport SEXP _DNAprofiles_ZcompKIpairswithtable(SEXP XSEXP, SEXP db1SEXP, SEXP db2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type db1(db1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type db2(db2SEXP);
    rcpp_result_gen = Rcpp::wrap(ZcompKIpairswithtable(X, db1, db2));
    return rcpp_result_gen;
END_RCPP
}
// ZcompKItargetsdbwithtable
NumericMatrix ZcompKItargetsdbwithtable(List X, IntegerMatrix db1, IntegerMatrix db2);
RcppExport SEXP _DNAprofiles_ZcompKItargetsdbwithtable(SEXP XSEXP, SEXP db1SEXP, SEXP db2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type db1(db1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type db2(db2SEXP);
    rcpp_result_gen = Rcpp::wrap(ZcompKItargetsdbwithtable(X, db1, db2));
    return rcpp_result_gen;
END_RCPP
}
// ZcompKIwithtable
NumericVector ZcompKIwithtable(List X, IntegerMatrix db);
RcppExport SEXP _DNAprofiles_ZcompKIwithtable(SEXP XSEXP, SEXP dbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type db(dbSEXP);
    rcpp_result_gen = Rcpp::wrap(ZcompKIwithtable(X, db));
    return rcpp_result_gen;
END_RCPP
}
// Zrmp
NumericMatrix Zrmp(IntegerMatrix db, NumericMatrix fr, double f, int retpermarker);
RcppExport SEXP _DNAprofiles_Zrmp(SEXP dbSEXP, SEXP frSEXP, SEXP fSEXP, SEXP retpermarkerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type db(dbSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type fr(frSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    Rcpp::traits::input_parameter< int >::type retpermarker(retpermarkerSEXP);
    rcpp_result_gen = Rcpp::wrap(Zrmp(db, fr, f, retpermarker));
    return rcpp_result_gen;
END_RCPP
}
// Zstl_nth_element
NumericVector Zstl_nth_element(NumericVector x, int n);
RcppExport SEXP _DNAprofiles_Zstl_nth_element(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(Zstl_nth_element(x, n));
    return rcpp_result_gen;
END_RCPP
}
// Zsumprodxy
double Zsumprodxy(NumericVector x, NumericVector y);
RcppExport SEXP _DNAprofiles_Zsumprodxy(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(Zsumprodxy(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DNAprofiles_Zprnextallele", (DL_FUNC) &_DNAprofiles_Zprnextallele, 4},
    {"_DNAprofiles_Zprnextalleles", (DL_FUNC) &_DNAprofiles_Zprnextalleles, 4},
    {"_DNAprofiles_Zallfinitepos", (DL_FUNC) &_DNAprofiles_Zallfinitepos, 3},
    {"_DNAprofiles_Zcountallelescpp", (DL_FUNC) &_DNAprofiles_Zcountallelescpp, 3},
    {"_DNAprofiles_Zdbcomparepairwise", (DL_FUNC) &_DNAprofiles_Zdbcomparepairwise, 3},
    {"_DNAprofiles_Zdbcomparepairwisetrackhits", (DL_FUNC) &_DNAprofiles_Zdbcomparepairwisetrackhits, 4},
    {"_DNAprofiles_Zdbcomparepairwisemc", (DL_FUNC) &_DNAprofiles_Zdbcomparepairwisemc, 4},
    {"_DNAprofiles_Zdbcomparepairwisemctrackhits", (DL_FUNC) &_DNAprofiles_Zdbcomparepairwisemctrackhits, 5},
    {"_DNAprofiles_Zdistapprox", (DL_FUNC) &_DNAprofiles_Zdistapprox, 5},
    {"_DNAprofiles_Zproductdist", (DL_FUNC) &_DNAprofiles_Zproductdist, 8},
    {"_DNAprofiles_Zexactq", (DL_FUNC) &_DNAprofiles_Zexactq, 6},
    {"_DNAprofiles_ZfindIntervalcpp", (DL_FUNC) &_DNAprofiles_ZfindIntervalcpp, 2},
    {"_DNAprofiles_Zki", (DL_FUNC) &_DNAprofiles_Zki, 11},
    {"_DNAprofiles_ZcompKIpairswithtable", (DL_FUNC) &_DNAprofiles_ZcompKIpairswithtable, 3},
    {"_DNAprofiles_ZcompKItargetsdbwithtable", (DL_FUNC) &_DNAprofiles_ZcompKItargetsdbwithtable, 3},
    {"_DNAprofiles_ZcompKIwithtable", (DL_FUNC) &_DNAprofiles_ZcompKIwithtable, 2},
    {"_DNAprofiles_Zrmp", (DL_FUNC) &_DNAprofiles_Zrmp, 4},
    {"_DNAprofiles_Zstl_nth_element", (DL_FUNC) &_DNAprofiles_Zstl_nth_element, 2},
    {"_DNAprofiles_Zsumprodxy", (DL_FUNC) &_DNAprofiles_Zsumprodxy, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_DNAprofiles(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
