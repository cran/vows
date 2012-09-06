\name{R2nii}
\alias{R2nii}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Save data to a NIfTI file
}
\description{
This function can be used to output the results of voxelwise RLRT or smoothing.}
\usage{
R2nii(arr, name.nii)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr}{a 3D or 4D array containing data to be saved.}
  \item{name.nii}{filename, excluding the .nii extension.}
}
\value{
None; a NIfTI file is created.}
\author{
Lei Huang \email{huangracer@gmail.com} and
Philip Reiss \email{phil.reiss@nyumc.org}
}

\seealso{
\code{\link{nii2R}}
}
