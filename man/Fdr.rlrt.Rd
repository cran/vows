\name{Fdr.rlrt}
\alias{Fdr.rlrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
False discovery rate estimation for massively parallel 
restricted likelihood ratio tests}
\description{
Given a set of RLRT results and a threshold, this function outputs an estimate of the FDR (in the empirical Bayes sense of Efron, 2010) when the given threshold is used to determine which null hypotheses to reject.}
\usage{
Fdr.rlrt(rlrt.obj, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rlrt.obj}{
an RLRT object obtained from \code{\link{rlrt.mp}} or \code{\link{rlrt4d}}.
}
  \item{threshold}{
threshold at which the null hypothesis is rejected.}
}

\value{A list with elements
\item{MoM}{FDR based on method of moments estimator of RLRT parameters (Greven et al., 2008).}
\item{ML}{FDR based on maximum likelihood estimation of RLRT parameters, as described in Greven et al. (2008).}
}
\references{
Efron, B. (2010).  \emph{Large-Scale Inference: Empirical Bayes Methods for Estimation, Testing, and Prediction}.  New York: Cambridge University Press.

Greven, S., Crainiceanu, C. M., Kuechenhoff, H., and Peters, A. (2008).  Restricted likelihood ratio testing for zero variance components in linear mixed models.  \emph{Journal of Computational and Graphical Statistics}, 17(4), 870--891.}
\author{
Philip Reiss \email{phil.reiss@nyumc.org}}

\seealso{
\code{\link{rlrt.mp}}, \code{\link{rlrt4d}}
}
\examples{
# See example for rlrt.mp
}

