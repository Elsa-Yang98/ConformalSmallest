% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridge_funs.R
\name{ginverse.fun}
\alias{ginverse.fun}
\title{Conformal prediction for linear regression}
\usage{
ginverse.fun(x, y, x0, alpha = 0.1)
}
\arguments{
\item{x}{A N*d training matrix}

\item{y}{A N*1 training vector}

\item{x0}{A N0*d testing vector}

\item{alpha}{miscoverage level}
}
\value{
upper and lower prediction intervals for X0
}
\description{
Conformal prediction for linear regression
}
