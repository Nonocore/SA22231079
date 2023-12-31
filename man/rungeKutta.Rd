% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rungeKutta}
\alias{rungeKutta}
\title{the fourth-order Runge-Kutta}
\usage{
rungeKutta(t0, y0, h, n, g, c, m)
}
\arguments{
\item{t0}{initial t value}

\item{y0}{initial y values}

\item{h}{step size}

\item{n}{number of iterations}

\item{g}{gravity acceleration}

\item{c}{air resistance coefficient}

\item{m}{object mass}
}
\value{
the result of the solution
}
\description{
the fourth-order Runge-Kutta method for solving a system of differential equations
}
\examples{
\dontrun{
rnC = rungeKutta(t0 = 0, y0 = 0, h = 0.0002, n = 10000, g = 9.8, c = 0.5, m = 1.0)
}
}
