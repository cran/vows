\name{lm.mp}
\alias{lm.mp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Massively parallel linear regression models
}
\description{
Efficiently fits \eqn{V} linear models with a common design matrix, where \eqn{V} may be very large, e.g., the number of voxels in a brain imaging application.}
\usage{
lm.mp(Y, formula, store.fitted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{\eqn{n \times V} outcome matrix.}
  \item{formula}{a formula object such as "\code{~ x1 + x2}".}
  \item{store.fitted}{logical: Should the fitted values be stored?  For large \eqn{V}, setting this to \code{TRUE} may cause memory problems.}
}

\value{
\item{coef}{\eqn{p \times V} matrix of coefficient estimates.}
\item{sigma2}{\eqn{V}-dimensional vector of error variance estimates.}
\item{se.coef}{\eqn{p \times V} matrix of coefficient standard error estimates.}
\item{X}{\eqn{n \times p} common design matrix.}
\item{fitted}{\eqn{n \times V} matrix of fitted values.}
}
\author{
Philip Reiss \email{phil.reiss@nyumc.org},
Lei Huang \email{huangracer@gmail.com}, and
Yin-Hsiu Chen \email{enjoychen0701@gmail.com}
}

\seealso{
\code{\link{lm4d}}, \code{\link{summary.lm.mp}}
}
\examples{
# Please see example for lm4d
}

