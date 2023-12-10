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

