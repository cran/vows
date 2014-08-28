#' Save data to a NIfTI file
#' 
#' This function can be used to output the results of voxelwise RLRT or
#' smoothing.
#' 
#' 
#' @param arr a 3D or 4D array containing data to be saved.
#' @param name.nii filename, excluding the .nii extension.
#' @return None; a NIfTI file is created.
#' @author Lei Huang \email{huangracer@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @seealso \code{\link{nii2R}}
#' @export
R2nii <-
function(arr, name.nii) {
	ndim = length(dim(arr))
    nim = Rniftilib::nifti.image.new()
    if (ndim==3) {
    	nim$dim = attributes(arr)$dim.nii
        arr[is.na(arr)] = 0
        nim[attributes(arr)$x.ind, attributes(arr)$y.ind, attributes(arr)$z.ind] = arr
    }
    else if (ndim==4) {
    	nim$dim = c(attributes(arr)$dim.nii, dim(arr)[4])
        arr[is.na(arr)] = 0
        nim[attributes(arr)$x.ind, attributes(arr)$y.ind, attributes(arr)$z.ind, ] = arr
    }
    Rniftilib::nifti.set.filenames(nim, name.nii)
    Rniftilib::nifti.image.write(nim)
}

