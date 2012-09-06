\name{plot.semipar.mp}
\alias{plot.semipar.mp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot massively parallel semiparametric models
}
\description{
Given a massively parallel smoothing object produced by \code{\link{semipar.mp}}, the function plots the fitted smooth(s) for a given point (e.g., at a given voxel).
}
\usage{
\method{plot}{semipar.mp}(x, Y, arr.ind = NULL, which.vox = NULL, which.smooth = NULL,
     coverage = 0.95, length.new = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "\code{\link{semipar.mp}}".
}
  \item{Y}{
an \eqn{n \times V} outcome matrix.
}
  \item{arr.ind}{
a 3-element vector specifying the element of the 3-dimensional array of locations (e.g., voxels) for which plotting is desired. If \code{NULL}, \code{which.vox} must be specified.
}
  \item{which.vox}{
the index of the voxel to be plotted. If \code{NULL}, \code{arr.ind} must be specified.
}
  \item{which.smooth}{
the index of the smooth term of which the confidence interval plot is to be displayed. The default value is \code{NULL} which refers to displaying the plots for all the smooth terms in the model. 
}
  \item{coverage}{
the confidence level of the pointwise confidence intervals in the plot.
}
  \item{length.new}{
length of the vector of ordered variables with which to predict.
}
  \item{\dots}{
arguments to be passed \code{\link[graphics]{plot}}
}
}

\author{
Yin-Hsiu Chen \email{enjoychen0701@gmail.com} and
Philip Reiss \email{phil.reiss@nyumc.org}
}



\examples{
n<-32
Ys <- matrix(0, n, 5)
for(i in 1:n) Ys[i,]<--2:2+rnorm(5, i^2, i^0.5)+sin(i)
x1 <- rnorm(n,0,5)
x2 <- 1:n+runif(n, 1, 20)
semipar.obj <- semipar.mp(~x1+sf(x2,k=10),Y=Ys,lsp=seq(5,50,,30))
plot(semipar.obj, Y=Ys, which.vox=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

