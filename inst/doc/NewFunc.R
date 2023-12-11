## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
#R：
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
f1 <- function(t, v) 9.8 - v^2
v1 <- RK4(dynamics = f1, y0 = 0, tn = 1, n = 5000)
plot(v1)

## ----rungeKutta()-------------------------------------------------------------
library(Rcpp)
ans <- rungeKutta(t0 = 0, y0 = c(0,0), h = 0.0002, n = 5000, g = 9.8, c = 1.0, m = 1.0)
plot(ans[,2])

## -----------------------------------------------------------------------------
library(microbenchmark)
tm1 <- microbenchmark::microbenchmark(
  rnR = RK4(dynamics = function(t, v) 9.8 - v^2, y0=0, t0 = 0, tn = 1, n = 5000),
  rnC = rungeKutta(t0 = 0, y0 = c(0,0), h = 0.0002, n = 5000, g = 9.8, c = 1.0, m = 1.0)
)
print(summary(tm1)[,c(1,3,5,6)])

## -----------------------------------------------------------------------------
#R：
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
f1 <- function(t, v) 9.8 - v^2
v1 <- EulerScheme(dynamics = f1, y0 = 0, tn = 1, n = 5000)
plot(v1)


## ----eulerMethod()------------------------------------------------------------
library(Rcpp)
ans <- eulerMethod(t0 = 0, y0 = c(0,0), h = 0.0002, n = 5000, g = 9.8, c = 1.0, m = 1.0)
plot(ans[,2])

## ----EulerScheme(),eulerMethod()----------------------------------------------
library(microbenchmark)
tm1 <- microbenchmark::microbenchmark(
  rnR = EulerScheme(dynamics = function(t, v) 9.8 - v^2, y0=0, t0 = 0, tn = 1, n = 5000),
  rnC = eulerMethod(t0 = 0, y0 = c(0,0), h = 0.0002, n = 5000, g = 9.8, c = 1, m = 1.0)
)
print(summary(tm1)[,c(1,3,5,6)])

