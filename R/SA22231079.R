#' @title Solve for first order ODEs by R-K 4th order
#' @description {One dimensional first order ODE solver using either RK4}
#' @param dynamics the function of t and y in the explicit ODE \eqn{y'=f(t,y)}
#' @param y0 the initial value
#' @param t0 the initial value input
#' @param tn the terminal value input
#' @param n the number of sub-intervals to use, defaults to 1000
#' @import ggplot2
#' @return data.frame of the time grid and solution grid
#' @examples
#' \dontrun{
#'   f1 <- function(t, v) 9.8 - v^2
#'   f2 <- function(t, v) 9.8 - 2*v^2
#'   f3 <- function(t, v) 9.8 - 3*v^2
#'   v1 <- RK4(dynamics = f1, y0 = 0, tn = 1, n = 5000)
#'   v2 <- RK4(dynamics = f2, y0 = 0, tn = 1, n = 5000)
#'   v3 <- RK4(dynamics = f3, y0 = 0, tn = 1, n = 5000)
#'   plot(v1, col = 'red')
#'   par(new=TRUE)
#'   plot(v2, col = 'blue')
#'   par(new=TRUE)
#'   plot(v3, col = 'green')
#' }
#' @export
RK4 <- function(dynamics, y0, t0 = 0, tn = 1, n = 1000)
{
  h <- (tn-t0)/n
  y <- matrix(0, n+1)
  # n+1 points, n sub-intervals
  tt <- seq(t0, tn, length.out = n+1)
  y[1] <- y0
  for(i in 2:(n+1))
  {
    k1 <- dynamics(tt[i-1], y[i-1])
    k2 <- dynamics(tt[i-1]+0.5*h, y[i-1]+0.5*h*k1)
    k3 <- dynamics(tt[i-1]+0.5*h, y[i-1]+0.5*h*k2)
    k4 <- dynamics(tt[i-1]+h, y[i-1]+h*k3)
    
    y[i] <- y[i-1] + (h/6)*(k1+2*k2+2*k3+k4)
  }
  return(data.frame(time = tt, state = y))
}


#' @title Solve for first order ODEs by EulerScheme
#' @description {One dimensional first order ODE solver using either EulerScheme}
#' @param dynamics the function of t and y in the explicit ODE \eqn{y'=f(t,y)}
#' @param y0 the initial value
#' @param t0 the initial value input
#' @param tn the terminal value input
#' @param n the number of sub-intervals to use, defaults to 1000
#' @import ggplot2
#' @return data.frame of the time grid and solution grid
#' @examples
#' \dontrun{
#'   f1 <- function(t, v) 9.8 - v^2
#'   f2 <- function(t, v) 9.8 - 2*v^2
#'   f3 <- function(t, v) 9.8 - 3*v^2
#'   v1 <- EulerScheme(dynamics = f1, y0 = 0, tn = 1, n = 5000)
#'   v2 <- EulerScheme(dynamics = f2, y0 = 0, tn = 1, n = 5000)
#'   v3 <- EulerScheme(dynamics = f3, y0 = 0, tn = 1, n = 5000)
#'   plot(v1, col = 'red')
#'   par(new=TRUE)
#'   plot(v2, col = 'blue')
#'   par(new=TRUE)
#'   plot(v3, col = 'green')
#' }
#' @export
EulerScheme <- function(dynamics, y0, t0 = 0, tn = 1, n = 1000)
{
  h <- (tn-t0)/n
  y <- matrix(0, n+1)
  # n+1 points, n sub-intervals
  tt <- seq(t0, tn, length.out = n+1)
  y[1] <- y0
  for(i in 2:(n+1))
  {
    y[i] <- y[i-1] + h*dynamics(tt[i-1], y[i-1])
  }
  return(data.frame(time = tt, state = y))
}

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C++ functions between RK4
#' @examples
#' \dontrun{
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = RK4(dynamics = function(t, v) 9.8 - v^2, y0=0, t0 = 0, tn = 1, n = 5000),
#'   rnC = rungeKutta(t0 = 0, y0 = 0, h = 0.0002, n = 10000, g = 9.8, c = 0.5, m = 1.0)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @useDynLib SA22231079
NULL

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C++ functions between Euler
#' @examples
#' \dontrun{
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = EulerScheme(dynamics = function(t, v) 9.8 - v^2, y0=0, t0 = 0, tn = 1, n = 5000),
#'   rnC = eulerMethod(t0 = 0, y0 = 0, h = 0.0002, n = 10000, g = 9.8, c = 0.5, m = 1.0)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @useDynLib SA22231079
NULL
