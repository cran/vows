\name{lm4d}
\alias{lm4d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Voxelwise linear models}
\description{
This is a wrapper function for \code{\link{lm.mp}} to handle 3D image responses.
}
\usage{
lm4d(arr4d, formula, store.fitted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr4d}{a 4-dimensional response array, where the first 3 dimensions refer to spatial coordinates and the last dimension corresponds to different images.}
  \item{formula,store.fitted}{see \code{\link{lm.mp}}.}
}

\value{
An object of class \code{"\link{lm.mp}"}, with two changes. (1) If \code{store.fitted = TRUE}, the fitted values are given as a 4-dimensional array. (2) A \code{call} component is included.
}
\author{
Lei Huang \email{huangracer@gmail.com},
Yin-Hsiu Chen \email{enjoychen0701@gmail.com}, and
Philip Reiss \email{phil.reiss@nyumc.org}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lm.mp}}
}
\examples{
data(test)
d4 = test$d4
x = test$x
lmobj = lm4d(d4, ~x)

# Convert d4 to a matrix, and confirm that lm.mp() gives the same results as lm4d()
d4.2 = d4
dim(d4.2) = c(prod(dim(d4)[1:3]), dim(d4)[4])
Y = t(d4.2)
lmobj2 = lm.mp(Y, ~x)
all.equal(lmobj$coef, lmobj2$coef)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
