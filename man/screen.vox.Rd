\name{screen.vox}
\alias{screen.vox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Screen voxels for a voxelwise smoothing object
}
\description{
Inputs a voxelwise smoothing object as produced by \code{\link{semipar4d}}, and outputs an object containing the results for a subset of the voxels.}
\usage{
screen.vox(semi.obj, arr4d, include)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{semi.obj}{an object of class \code{\link{semipar.mp}}.}
  \item{arr4d}{the 4-dimensional array used to generate the object.}
  \item{include}{a logical matrix indicating which points (or voxels) should be included.}
}

\value{
a modified version of \code{semipar.obj}, with pointwise coefficients (\code{coef} component), pointwise degrees of freedom (\code{pwdf}), pointwise log smoothing parameter (\code{pwlsp}), and pointwise variance estimate (\code{sigma2}) for the points specified by \code{include} only.  
}
\author{
Lei Huang \email{huangracer@gmail.com} and
Philip Reiss \email{phil.reiss@nyumc.org}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{semipar.mp}}
}
\examples{
data(test)
d4 = test$d4
x = test$x
vw.obj = semipar4d(d4, formula = ~sf(x), data = data.frame(x = x), lsp=-5:5)

# Include only the first 600 voxels
sv = screen.vox(vw.obj, d4, rep(1:0, c(600,400)))}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
