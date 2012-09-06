\name{funkpanel}
\alias{funkpanel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive visualization of massively parallel smoothing results 
}
\description{
This function uses interactive graphics tools, provided by the \pkg{rpanel} package (Bowman et al., 2007), to visualize the results of functional k-means clustering as implemented by \code{\link{funkmeans}}.}
\usage{
funkpanel(fkmobj, semiobj, arr4d, predictor, titl = "", xlab = "",
          ylab = "", ncluster = nrow(fkmobj$centers),
          slice = dim(fkmobj$arr.cluster)[3]\%/\%2, ylim.scatter = NULL,
          deriv.legend = 0, ylim.legend = NULL, scattermain = NULL,
          colvec = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fkmobj}{
a functional k-means clustering object obtained from \code{\link{funkmeans4d}}.
}
  \item{semiobj}{
the massively parallel smoothing object on which the clustering was based; ordinarily produced by \code{\link{semipar.mp}} or \code{\link{semipar4d}}.
}
  \item{arr4d}{
a 4-dimensional array containing the raw data that were smoothed at each point.  The first 3 dimensions refer to x, y, and z coordinates and the last dimension corresponds to different images.}
  \item{predictor}{
a vector or matrix of covariates.
}
  \item{titl}{
title of the panel.
}
  \item{xlab,ylab}{
x- and y-axis labels.
}
  \item{ncluster}{
number of clusters to display.  By default, all are displayed.
}
  \item{slice}{
index of the slice to be shown initially in the panel.
}
  \item{ylim.scatter}{
the y limits of the voxelwise scatterplots.
}
  \item{deriv.legend}{
which derivative to plot in the "legend"; see Details.  By default, the curves themselves are used.
}
  \item{ylim.legend}{
the y limits used in the "legend"; see Details.
}
  \item{scattermain}{
title for the scatter plots.
}
  \item{colvec}{
a vector of colors for the clusters. By default, this is set to the first \code{ncluster} elements of \code{c("dodgerblue", "green", "red", "orange", "yellow", "orchid",}  \code{ "brown", "grey", "purple")}, if \code{ncluster < 9}. 
}
}
\details{
The clusters defined by the input object are shown for a cross-section (slice), and a slider allows viewing of different slices.  Clicking on a voxel produces a scatterplot of the data at that voxel, along with the fitted scatterplot.  The "legend", produced by \code{\link{plot.funkmeans}}, consists of plots of 30 randomly selected curves, along with the cluster center, from each cluster.}
\references{
	Bowman, A., Crawford, E., Alexander, G., and Bowman, R. (2007). rpanel: Simple interactive controls for R functions using the tcltk package. \emph{Journal of Statistical Software}, 17(9).
}
\author{
Lei Huang \email{huangracer@gmail.com},
Yin-Hsiu Chen \email{enjoychen0701@gmail.com}, and Lan Huo \email{lan.huo@nyumc.org}
}

\seealso{
\code{\link{funkmeans}}, \code{\link{funkmeans4d}}, \code{\link{plot.funkmeans}}
}
\examples{
data(test)
d4 = test$d4
x = test$x
semi.obj = semipar4d(d4, ~sf(x), -5:5, data.frame(x = x))
fkmobj = funkmeans4d(semi.obj, d4, ncomp=6, centers=3)
funkpanel(fkmobj, semi.obj, d4, x)
}

