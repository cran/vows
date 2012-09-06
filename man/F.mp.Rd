\name{F.mp}
\alias{F.mp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
F-tests for massively parallel linear models}
\description{
Performs F-tests for removing one or more terms from each of a large number of models with common design matrix.}
\usage{
F.mp(formula, which)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula such as "\code{Y ~ X}", where \code{Y} is an \eqn{n \times V} response matrix and \code{X} is an \eqn{n \times p} design matrix common to all \eqn{V} models.
}
  \item{which}{
number or vector indicating which column(s) of the model matrix are to be tested for removal from the model.}
}
\value{
\item{F}{F-statistics for each of the models.}
\item{df1}{numerator degrees of freedom.}
\item{df2}{denominator degrees of freedom.}
\item{pvalue}{upper-tailed p-value.}
\item{X}{design matrix.}
}

\author{
Philip Reiss \email{phil.reiss@nyumc.org} and Lei Huang \email{huangracer@gmail.com}
}

\seealso{
\code{\link{lm.mp}}, \code{\link{permF.mp}}
}
\examples{
Y = matrix(rnorm(6000), nrow=20)
X = rnorm(20)
t2 = F.mp(Y~X, which=2)
}