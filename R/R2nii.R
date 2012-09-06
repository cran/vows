R2nii <-
function(arr, name.nii) {
	ndim = length(dim(arr))
    require(Rniftilib)
    nim = nifti.image.new()
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
    nifti.set.filenames(nim, name.nii)
    nifti.image.write(nim)
}

