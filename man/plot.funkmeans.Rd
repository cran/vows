\name{plot.funkmeans}
\alias{plot.funkmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting of k-means clustering results for massively parallel smooths
}
\description{
Visualization of functional k-means clustering as implemented by \code{\link{funkmeans}}.}
\usage{
\method{plot}{funkmeans}(x, fdobj, deriv = 0, ylim = NULL, ncluster = nrow(x$centers), mfrow = NULL,
     colvec = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a functional k-means clustering object obtained from \code{\link{funkmeans}}.
}
  \item{fdobj}{
a functional data object, of class \code{"\link[fda]{fd}"}, defining the set of curves being clustered. See the example.
}
  \item{deriv}{
which derivative to display in the plots, which show 30 randomly selected curves, along with the cluster center, from each cluster.  By default, the "0th derivative" is used (i.e., the curves themselves).
}
  \item{ylim}{
the y limits for the plots.
}
  \item{ncluster}{
number of clusters to display.  By default, all are displayed.
}
\item{mfrow}{a vector of length 2 giving the numbers of rows and columns for the array of plots. By default, the number of rows will exceed the number of columns by \code{0} or \code{1}, depending on \code{ncluster}.}
  \item{colvec}{
a vector of colors for the clusters. By default, this is set to the first \code{ncluster} elements of \code{c("dodgerblue", "green", "red", "orange", "yellow", "orchid",}  \code{ "brown", "grey", "purple")}, if \code{ncluster <= 9}. 
}
  \item{\dots}{
arguments passed to \code{\link[graphics]{plot}.}
}
}

\author{
Yin-Hsiu Chen \email{enjoychen0701@gmail.com},
Philip Reiss \email{phil.reiss@nyumc.org}
and Lan Huo \email{lan.huo@nyumc.org}
}

\seealso{
\code{\link{funkmeans}}, \code{\link{funkpanel}}
}


\examples{
data(test)
d4 = test$d4
x = test$x
semi.obj = semipar4d(d4, formula = ~sf(x), data = data.frame(x = x), lsp=-5:5)
fdobj = fd(coef = semi.obj$coef[-1, ], basis = semi.obj$list.all[[1]]$basis)
fkmobj = funkmeans(fdobj, ncomp = 8, centers = 6)
plot(fkmobj, fdobj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

