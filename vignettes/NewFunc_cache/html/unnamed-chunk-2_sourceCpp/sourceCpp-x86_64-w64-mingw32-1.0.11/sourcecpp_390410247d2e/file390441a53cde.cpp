#include <Rcpp.h>
using namespace Rcpp;

// Define the system of differential equations
void f(double t, NumericVector y, NumericVector dydt, double g, double c, double m) {
  dydt[0] = y[1];
  dydt[1] = g - (c / m) * y[1] * y[1];
}

// [[Rcpp::export]]
NumericMatrix rungeKutta(double t0, NumericVector y0, double h, int n, double g, double c, double m) {
  int dim = y0.length();
  
  NumericMatrix y(n + 1, dim);
  NumericVector k1(dim), k2(dim), k3(dim), k4(dim);
  
  double t = t0;
  y(0, _) = y0;
  
  for (int i = 0; i < n; i++) {
    f(t, y(i, _), k1, g, c, m);
    f(t + h / 2, y(i, _) + h * k1 / 2, k2, g, c, m);
    f(t + h / 2, y(i, _) + h * k2 / 2, k3, g, c, m);
    f(t + h, y(i, _) + h * k3, k4, g, c, m);
    
    y(i + 1, _) = y(i, _) + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    
    t += h;
  }
  
  return y;
}



#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rungeKutta
NumericMatrix rungeKutta(double t0, NumericVector y0, double h, int n, double g, double c, double m);
RcppExport SEXP sourceCpp_1_rungeKutta(SEXP t0SEXP, SEXP y0SEXP, SEXP hSEXP, SEXP nSEXP, SEXP gSEXP, SEXP cSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(rungeKutta(t0, y0, h, n, g, c, m));
    return rcpp_result_gen;
END_RCPP
}
