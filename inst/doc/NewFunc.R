## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  #R：
#  RK4 <- function(dynamics, y0, t0 = 0, tn = 1, n = 1000)
#  {
#    h <- (tn-t0)/n
#    y <- matrix(0, n+1)
#    # n+1 points, n sub-intervals
#    tt <- seq(t0, tn, length.out = n+1)
#    y[1] <- y0
#    for(i in 2:(n+1))
#    {
#      k1 <- dynamics(tt[i-1], y[i-1])
#      k2 <- dynamics(tt[i-1]+0.5*h, y[i-1]+0.5*h*k1)
#      k3 <- dynamics(tt[i-1]+0.5*h, y[i-1]+0.5*h*k2)
#      k4 <- dynamics(tt[i-1]+h, y[i-1]+h*k3)
#  
#      y[i] <- y[i-1] + (h/6)*(k1+2*k2+2*k3+k4)
#    }
#    return(data.frame(time = tt, state = y))
#  }
#  f1 <- function(t, v) 9.8 - v^2
#  v1 <- RK4(dynamics = f1, y0 = 0, tn = 1, n = 5000)
#  #C++
#  # library(Rcpp)
#  # #include <Rcpp.h>
#  # using namespace Rcpp;
#  #
#  # // [[Rcpp::export]]
#  # void f(double t, NumericVector y, NumericVector dydt, double g, double c, double m) {
#  #   dydt[0] = y[1];
#  #   dydt[1] = g - (c / m) * y[1] * y[1];
#  # }
#  # // [[Rcpp::export]]
#  # NumericMatrix rungeKutta(double t0, NumericVector y0, double h, int n, double g, double c, double m) {
#  #   int dim = y0.length();
#  #
#  #   NumericMatrix y(n + 1, dim);
#  #   NumericVector k1(dim), k2(dim), k3(dim), k4(dim);
#  #
#  #   double t = t0;
#  #   y(0, _) = y0;
#  #
#  #   for (int i = 0; i < n; i++) {
#  #     f(t, y(i, _), k1, g, c, m);
#  #     f(t + h / 2, y(i, _) + h * k1 / 2, k2, g, c, m);
#  #     f(t + h / 2, y(i, _) + h * k2 / 2, k3, g, c, m);
#  #     f(t + h, y(i, _) + h * k3, k4, g, c, m);
#  #
#  #     y(i + 1, _) = y(i, _) + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
#  #
#  #     t += h;
#  #   }
#  #
#  #   return y;
#  # }
#  

## ----eval=FALSE---------------------------------------------------------------
#  #R：
#  EulerScheme <- function(dynamics, y0, t0 = 0, tn = 1, n = 1000)
#  {
#    h <- (tn-t0)/n
#    y <- matrix(0, n+1)
#    # n+1 points, n sub-intervals
#    tt <- seq(t0, tn, length.out = n+1)
#    y[1] <- y0
#    for(i in 2:(n+1))
#    {
#      y[i] <- y[i-1] + h*dynamics(tt[i-1], y[i-1])
#    }
#    return(data.frame(time = tt, state = y))
#  }
#  f1 <- function(t, v) 9.8 - v^2
#  v1 <- EulerScheme(dynamics = f1, y0 = 0, tn = 1, n = 5000)
#  # C++
#  # library(Rcpp)
#  # #include <Rcpp.h>
#  # using namespace Rcpp;
#  # // [[Rcpp::export]]
#  # void f2(double t, NumericVector y, NumericVector dydt, double g, double c, double m) {
#  #   dydt[0] = y[1];
#  #   dydt[1] = g - (c / m) * y[1] * y[1];
#  # }
#  # // [[Rcpp::export]]
#  # NumericMatrix eulerMethod(double t0, NumericVector y0, double h, int n, double g, double c, double m) {
#  #   int dim = y0.length();
#  #
#  #   NumericMatrix y(n + 1, dim);
#  #
#  #   double t = t0;
#  #   y(0, _) = y0;
#  #
#  #   for (int i = 0; i < n; i++) {
#  #     NumericVector dydt(dim);
#  #     f2(t, y(i, _), dydt, g, c, m);
#  #
#  #     y(i + 1, _) = y(i, _) + h * dydt;
#  #
#  #     t += h;
#  #   }
#  #
#  #   return y;
#  # }
#  
#  

