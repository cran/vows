rlr.xz <-
function(x, nbasis=15, norder=4, nulldim=NULL, B=NULL, P) {
    if (is.null(B))   {
        bsb = create.bspline.basis(range(x, na.rm=TRUE), nbasis, norder)
        B = eval.basis(x, bsb)
        if (is.null(nulldim)) nulldim = norder-2
        P = getbasispenalty(bsb, nulldim)
    }
    else nbasis = ncol(B)
    svd.pen = svd(P)
    Us = svd.pen$u[ , 1:(nbasis-nulldim)]; Un = svd.pen$u[ , (nbasis-nulldim+1):nbasis]  
    X = B %*% Un
    Z = scale(B %*% Us, FALSE, sqrt(svd.pen$d[1:(nbasis-nulldim)]))  # see Wood (2004)
    list(X=X, Z=Z)
}

