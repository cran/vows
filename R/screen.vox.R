screen.vox <-
function(semi.obj, arr4d, include)    {
    has.data = attributes(arr4d)$has.data
    data.inds = which(has.data==TRUE, arr.ind=TRUE)
    if (ncol(semi.obj$coef)!=length(include))  stop("Dimensions don't match!")
    semi.obj$coef = semi.obj$coef[ , include]
    semi.obj$pwdf = semi.obj$pwdf[include]
    semi.obj$pwlsp = semi.obj$pwlsp[include]
    semi.obj$sigma2 = semi.obj$sigma2[include]
    semi.obj$include = has.data
    semi.obj$include[has.data & !is.na(has.data)] = include
    return(semi.obj)
}

