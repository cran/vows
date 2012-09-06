\name{funkmeans}
\alias{funkmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional k-means clustering for parallel smooths}
\description{
This function performs k-means clustering for curve estimates corresponding to each of a 3D grid of points.
For example, when scatterplot smoothing is performed at each of a grid of brain voxels as in Reiss et al. (2013), this function can be used to cluster the obtained smooths.}
\usage{
funkmeans(obj, which.smooth = 1, deriv = 1, lambda = 0, ncomp, centers, nstart = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a "massively parallel" semiparametric model object, of class \code{semipar.mp}. 
}
  \item{which.smooth}{which of the smooth terms to use for functional k-means: e.g., if \code{obj} consists of fits of the model \code{~ sf(x1) + sf(x2)}, \code{which.smooth = 2} yields k-means clustering for the curves with respect to \code{x2}. The default suffices for models with a single smooth term.
}
  \item{deriv}{which derivative of the curves should be clustered.  If \code{0}, the curves themselves are clustered; if \code{1} (the default), their first derivatives are clustered, a natural way to assign curves of similar shape to the same cluster.
}
 \item{lambda}{smoothing parameter for functional PCA as implemented by \code{\link[fda]{pca.fd}}.
}
  \item{ncomp}{number of functional principal components.
}
  \item{centers}{number of clusters.
}
  \item{nstart}{number of randomly chosen sets of initial centers used by the \code{\link[stats]{kmeans}} function.
}
}
\details{
The functional clustering algorithm consists of performing (i) functional principal component analysis of the curve estimates or their derivatives, followed by (ii) k-means clustering of the functional PC scores (Tarpey and Kinateder, 2003).}
\value{
An object of class "funkmeans", which is a list with elements:
\item{cluster, centers, withinss, tots, tot.withinss, betweenness, size}{see \code{\link[stats]{kmeans}}.}
\item{basis,coef}{basis object and coefficient matrix defining the functional data object (see \code{\link[fda]{fd}}) for the curves that are clustered.}
\item{fpca}{functional principal components object, output by \code{\link[fda]{pca.fd}}.}
\item{R2}{proportion of variance explained by the k clusters.}
\item{which.smooth}{value of argument \code{which.smooth} in the function call.}
}
\references{
Reiss, P. T., Huang, L., Chen, Y.-H., Huo, L., Tarpey, T., and Mennes, M. (2013). Massively parallel nonparametric regression, with an application to developmental brain mapping. \emph{Journal of Computational and Graphical Statistics}, to appear.  Available at \url{http://works.bepress.com/phil_reiss/24}

Tarpey, T., and Kinateder, K. K. J. (2003).  Clustering functional data.  \emph{Journal of Classification}, 20, 93--114.
}
\author{
Philip Reiss \email{phil.reiss@nyumc.org}, Lei Huang \email{huangracer@gmail.com}
and Lan Huo \email{lan.huo@nyumc.org}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{funkmeans4d}}, \code{\link{funkpanel}}
}
\examples{
# See example for funkpanel
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
