\name{sf}
\alias{sf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Defining smooth functions in semiparametric model formulae
}
\description{
This function is called by \code{\link{semipar.mp}} to define B-spline smooths.
}
\usage{
sf(argvals, effect = NULL, k = 10, norder = 4, pen.order = 2, range.basis = NULL, 
   knots = "quantile")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{argvals}{
a vector or matrix of covariates.
}
  \item{effect}{
predictor whose effect varies with respect to \code{argvals}. E.g., if the effect of \code{diagnosis} varies with \code{age}, use \code{sf(age, effect = diagnosis)}. Similar to argument \code{by} in \code{\link[mgcv]{s}}.
}
  \item{k}{
number of B-spline basis functions.
}
  \item{norder}{
order of B-splines: the default, \code{4}, gives cubic B-splines.
}
  \item{pen.order}{
order of the penalty, i.e., of the derivative defining the penalty.
}
  \item{range.basis}{
a numeric vector of length 2 defining the interval over which the B-spline basis is created. If \code{NULL}, set to the range of the variable. 
}
  \item{knots}{
knots placement method for B-spline smoothing. The default, "quantile", places the knots at equally spaced quantiles of the data; "equispaced" gives equally spaced knots.}
}


\author{
Yin-Hsiu Chen \email{enjoychen0701@gmail.com} and
Philip Reiss \email{phil.reiss@nyumc.org}
}

