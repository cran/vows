% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funkmeans4d.R
\name{funkmeans4d}
\alias{funkmeans4d}
\title{Functional k-means clustering for parallel smooths for 4-dimensional data}
\usage{
funkmeans4d(fdobj, arr4d, ...)
}
\arguments{
\item{fdobj}{a functional data object, of class \code{"\link[fda]{fd}"},
defining the set of curves being clustered.}

\item{arr4d}{a 4-dimensional array containing the raw data that were
smoothed at each point.  The first 3 dimensions refer to x, y, and z
coordinates and the last dimension corresponds to different images.}

\item{\dots}{other arguments, passed to \code{\link{funkmeans}}.}
}
\value{
An object of class "funkmeans4d", which is also of class
\code{"\link{funkmeans}"} but has the additional component
\code{arr.cluster}: an array, of dimension \code{dim(arr4d)[1:3]}, giving
the cluster memberships.
}
\description{
This is a wrapper function for \code{\link{funkmeans}} to handle 3D image
responses.
}
\examples{

# See example for funkmeans
}
\author{
Philip Reiss \email{phil.reiss@nyumc.org}, Lei Huang
\email{huangracer@gmail.com} and Lan Huo
}
\seealso{
\code{\link{funkmeans}}
}

