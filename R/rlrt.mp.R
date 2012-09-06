rlrt.mp <-
function(Y, x=NULL, loginvsp, nbasis = 15, norder=4, nulldim=NULL, evalarg=NULL, get.df=FALSE, B=NULL, P=NULL) {
  	require(fda)
    if (is.null(x) & is.null(B))   stop("You must enter either x or B")
    if (is.null(x) & is.null(P))   stop("You must enter P")
    
    n = if (!is.null(x)) NROW(x) else NROW(B)
  	
    xz = rlr.xz(x, nbasis=nbasis, norder=norder, nulldim=nulldim, B=B, P=P)
    rlobj = rlrt.mp.fit(Y=Y, X=xz$X, Z=xz$Z, loginvsp=loginvsp, evalarg=evalarg, get.df=get.df)
    rlobj
}

