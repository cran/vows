\name{funkmeans4d}
\alias{funkmeans4d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional k-means clustering for parallel smooths for 4-dimensional data}
\description{
This is a wrapper function for \code{\link{funkmeans}} to handle 3D image responses.}
\usage{
funkmeans4d(fdobj, arr4d, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdobj}{a functional data object, of class \code{"\link[fda]{fd}"}, defining the set of curves being clustered.
}
  \item{arr4d}{a 4-dimensional array containing the raw data that were smoothed at each point.  The first 3 dimensions refer to x, y, and z coordinates and the last dimension corresponds to different images.
}
  \item{\dots}{
other arguments, passed to \code{\link{funkmeans}}.
}
}

\value{
An object of class "funkmeans4d", which is also of class \code{"\link{funkmeans}"} but has the additional component \code{arr.cluster}: an array, of dimension \code{dim(arr4d)[1:3]}, giving the cluster memberships.
}
\author{
Philip Reiss \email{phil.reiss@nyumc.org}, Lei Huang \email{huangracer@gmail.com}
and Lan Huo \email{lan.huo@nyumc.org}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{funkmeans}}, \code{\link{funkpanel}}
}
\examples{
# See example for funkpanel
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
