funkmeans4d <- function(obj, arr4d, ...) {
    has.data = attributes(arr4d)$has.data
    fkmobj <- funkmeans(obj=obj, ...)
  	# centers.fdobj = fd(coef=fkmobj$fpca$meanfd$coef %*% matrix(1,1,centers) + fkmobj$fpca$harmonics$coef %*% t(fkmobj$centers), basisobj=temp$basis)
  	include = if (is.null(obj$include)) has.data else obj$include
  	arr.cluster = array(NA, dim(include))
  	arr.cluster[include & !is.na(include)] = fkmobj$cluster
  	arr.cluster[has.data & !is.na(include) & !include] = 0
	attr(arr.cluster, "x.ind") = attr(arr4d, "x.ind")
	attr(arr.cluster, "y.ind") = attr(arr4d, "y.ind")
	attr(arr.cluster, "z.ind") = attr(arr4d, "z.ind")
	attr(arr.cluster, "dim.nii") = attr(arr4d, "dim.nii")
	# May need to get some other attributes from arr4d
    # fkmobj$centers.fdobj = centers.fdobj
    fkmobj$arr.cluster = arr.cluster
    class(fkmobj) = c(class(fkmobj), "funkmeans4d")
    fkmobj
}

