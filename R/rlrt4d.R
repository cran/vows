rlrt4d <-
function(arr4d, x=NULL, nbasis = 15, norder=4, nulldim=NULL, loginvsp, get.df=FALSE, B=NULL, P=NULL) {
    ndim = length(dim(arr4d))
    has.data <- attributes(arr4d)$has.data
    N = NROW(x)
    dim(arr4d) = c(prod(dim(arr4d)[1:(ndim-1)]), N)
    Y.d = t(arr4d[(as.vector(has.data)), ])
    rm(arr4d)  # i.e., remove from current environment
    rlobj <- rlrt.mp(Y.d, x=x, nbasis = nbasis, norder=norder, nulldim=nulldim, loginvsp=loginvsp, get.df=get.df, B=B, P=P)
    rlobj$call = match.call()
    class(rlobj) = "rlrt4d"
    return(rlobj)
}
