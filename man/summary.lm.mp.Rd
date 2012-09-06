\name{summary.lm.mp}
\alias{summary.lm.mp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing massively parallel linear model fits}
\description{
\code{summary} method for class "\code{lm.mp}".
}
\usage{
\method{summary}{lm.mp}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{lm.mp}}, ordinarily created by the function of that name or by \code{\link{lm4d}}.}
  \item{\dots}{not currently used.}
}

\value{
\item{tstat}{matrix of pointwise t-statistics for each coefficient in the linear model}
\item{pvalue}{matrix of the pointwise p-values for each coefficient in the linear model}
\item{aicc}{vector of pointwise corrected AIC}
}
\author{
Philip Reiss \email{phil.reiss@nyumc.org} and Lei Huang \email{huangracer@gmail.com}
}
\seealso{
\code{\link{lm.mp}}
}
\examples{
Y = matrix(rnorm(6000), nrow=20)
X = rnorm(20)
t1 = lm.mp(Y, ~X)
st1 = summary(t1)
}

