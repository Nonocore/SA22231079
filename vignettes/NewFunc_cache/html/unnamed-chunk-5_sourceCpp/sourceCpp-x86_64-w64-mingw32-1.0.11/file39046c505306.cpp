#include <Rcpp.h>
using namespace Rcpp;


// Define the system of differential equations
void f2(double t, NumericVector y, NumericVector dydt, double g, double c, double m) {
  dydt[0] = y[1];
  dydt[1] = g - (c / m) * y[1] * y[1];
}

// [[Rcpp::export]]
NumericMatrix eulerMethod(double t0, NumericVector y0, double h, int n, double g, double c, double m) {
  int dim = y0.length();
  
  NumericMatrix y(n + 1, dim);
  
  double t = t0;
  y(0, _) = y0;
  
  for (int i = 0; i < n; i++) {
    NumericVector dydt(dim);
    f2(t, y(i, _), dydt, g, c, m);
    
    y(i + 1, _) = y(i, _) + h * dydt;
    
    t += h;
  }
  
  return y;
}
