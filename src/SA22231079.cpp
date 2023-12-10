#include <Rcpp.h>
using namespace Rcpp;

// Define the system of differential equations
void f(double t, NumericVector y, NumericVector dydt, double g, double c, double m) {
  dydt[0] = y[1];
  dydt[1] = g - (c / m) * y[1] * y[1];
}

//' @title the fourth-order Runge-Kutta 
//' @description the fourth-order Runge-Kutta method for solving a system of differential equations
//' @param t0 initial t value
//' @param y0 initial y values
//' @param h step size
//' @param n number of iterations
//' @param g gravity acceleration
//' @param c air resistance coefficient
//' @param m object mass
//' @return the result of the solution
//' @examples
//' \dontrun{
//' rnC = rungeKutta(t0 = 0, y0 = 0, h = 0.0002, n = 10000, g = 9.8, c = 0.5, m = 1.0)
//' }
//' @export
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




// Define the system of differential equations
void f2(double t, NumericVector y, NumericVector dydt, double g, double c, double m) {
  dydt[0] = y[1];
  dydt[1] = g - (c / m) * y[1] * y[1];
}

//' @title the fourth-order eulerMethod
//' @description the fourth-order eulerMethod method for solving a system of differential equations
//' @param t0 initial t value
//' @param y0 initial y values
//' @param h step size
//' @param n number of iterations
//' @param g gravity acceleration
//' @param c air resistance coefficient
//' @param m object mass
//' @return the result of the solution
//' @examples
//' \dontrun{
//' rnC = eulerMethod(t0 = 0, y0 = 0, h = 0.0002, n = 10000, g = 9.8, c = 0.5, m = 1.0)
//' }
//' @export
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


