\name{rlrt4d}
\alias{rlrt4d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Voxelwise restricted likelihood ratio tests}
\description{
A wrapper function for \code{\link{rlrt.mp}} to handle 3D image responses.
}
\usage{
rlrt4d(arr4d, x = NULL, nbasis = 15, norder = 4, nulldim = NULL,
       loginvsp, get.df = FALSE, B = NULL, P = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr4d}{a 4-dimensional response array, where the first 3 dimensions refer to spatial coordinates and the last dimension corresponds to different images.}
  \item{x, nbasis, norder, nulldim, loginvsp, get.df, B, P}{
see \code{\link{rlrt.mp}}.
}
}
\value{
A massively parallel RLRT object, as produced by \code{\link{rlrt.mp}}.
}
\author{
Lei Huang \email{huangracer@gmail.com} and
Philip Reiss \email{phil.reiss@nyumc.org}
}


\seealso{
\code{\link{plot.rlrt4d}}, \code{\link{rlrt.mp}}
}
\examples{
data(test)
d4 = test$d4
x = test$x
rlrtobj = rlrt4d(d4, x, loginvsp = -5:5)
plot(rlrtobj, d4, slice=5)
\donttest{rlrtpanel(rlrtobj, d4, x)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
