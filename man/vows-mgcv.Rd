% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vows-mgcv.R
\name{vows-mgcv}
\alias{Predict.matrix.bspline.smooth}
\alias{smooth.construct.be.smooth.spec}
\alias{smooth.construct.bq.smooth.spec}
\alias{vows-mgcv}
\title{Utility functions related to the mgcv package}
\usage{
\method{smooth.construct}{bq.smooth.spec}(object, data, knots)

\method{smooth.construct}{be.smooth.spec}(object, data, knots)

\method{Predict.matrix}{bspline.smooth}(object, data)
}
\arguments{
\item{object}{a \code{\link[mgcv]{gam}} smooth specification object
generated by a term such as \code{s(x,bs="be")} or \code{s(x,bs="bq")}.}

\item{data}{For \code{smooth.construct.be.smooth.spec} and
\code{smooth.construct.bq.smooth.spec}, a list containing just the data
(including any \code{by} variable) required by the given term, with names
corresponding to \code{object$term} (and \code{object$by}). The \code{by}
variable is the last element. For \code{Predict.matrix.bspline.smooth}, a
data frame containing the values of the (named) covariates at which the
smooth term is to be evaluated. Exact requirements are as for
\code{\link[mgcv]{smooth.construct}} and
\code{\link[mgcv]{smooth.construct2}}.}

\item{knots}{a list containing any knots supplied for basis setup, in the
same order and with the same names as \code{data}. If \code{NULL}, a default
set of knots is used.}
}
\value{
Either \code{smooth.construct.be.smooth.spec} or
\code{smooth.construct.bq.smooth.spec} produces an object of class
\code{"bspline.smooth"}; see \code{\link[mgcv]{smooth.construct}} for the
elements that this object will contain. \code{Predict.matrix.bspline.smooth}
produces a matrix mapping the coefficients for the smooth term to its values
at the supplied data values.
}
\description{
These internal functions are used by \code{\link{semipar.mix.mp}} (but can
also be used more generally) to customize the implementation of B-spline
smoothing by \code{\link[mgcv]{gam}}. Specifically, a B-spline smooth with
equispaced knots can be incorporated in a call to \code{\link[mgcv]{gam}}
using a term of the form \code{s(x,bs="be")}, whereas knots at equally
spaced quantiles of the data can be specified by \code{s(x,bs="bq")}.
}
\details{
These functions are not normally called directly. For further details,
please see \code{\link[mgcv]{smooth.construct.ps.smooth.spec}} and
\code{\link[mgcv]{Predict.matrix.cr.smooth}}.
}
\examples{

x. = rnorm(20)
smoo.be <- smooth.construct.be.smooth.spec(s(x), data.frame(x = x.), NULL)
smoo.bq <- smooth.construct.bq.smooth.spec(s(x), data.frame(x = x.), NULL)
Predict.matrix.bspline.smooth(smoo.bq, data.frame(x = seq(min(x.),max(x.),,100)))

}
\author{
Yin-Hsiu Chen \email{enjoychen0701@gmail.com} and Philip Reiss
\email{phil.reiss@nyumc.org}
}

