sf <- function(argvals, effect=NULL, k = 10, norder = 4, pen.order = 2, range.basis = NULL, knots = "quantile") {
	require(fda)
    if (is.null(range.basis)) range.basis = range(argvals)
    if (knots == "quantile")	basis = create.bspline.basis(range.basis, breaks = quantile(argvals, seq(0,1,length.out=k-norder+2)), norder=norder)
    else  if (knots == "equispaced") basis = create.bspline.basis(range.basis, norder=norder, nbasis = k)
	modmat = eval.basis(argvals, basis) 
	if (!is.null(effect)) modmat = diag(effect) %*% modmat
    penmat = getbasispenalty(basis, pen.order)   
    constraint = if (is.null(effect)) colSums(modmat) else NULL  
    
    sf.out = list(basis = basis, modmat = modmat, penmat = penmat, 
                  constraint = constraint, argvals = argvals, effect = effect, k = k, 
                  norder = norder, pen.order = pen.order)
    class(sf.out) = "sf"
    sf.out
}

