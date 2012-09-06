# Must change pwdf; see end of Wood (2004), sec. 1.2. 
# Also, should add standard error estimates!  (PR, Dec. 29/11)
qplsc.mp <- function(Y, modmat, penmat, constr.list=NULL, lsp, nulldim=NULL, store.reml=FALSE, store.fitted=FALSE) {
   	n = nrow(Y)
   	if (nrow(modmat) != n) stop(paste("number of rows of 'x' or 'modmat' must match last dimension of 'Y'"))
  	L = NCOL(Y)
  	if (L==1) stop("'Y' must be a matrix with at least 2 columns.")
  	if (is.null(constr.list)) { 
        modmat.u <- modmat
        penmat.u <- penmat
        ttu <- diag(1, ncol(modmat))
    }  
    else {
        # Transform to unconstrained problem
        # ttu = qr.Q(qr(t(constrmat)), complete=TRUE)[ , -(1:NROW(constrmat))]
        ttu.ful <- diag(1, ncol(modmat))
        for (i in 1:length(constr.list)) 
          	ttu.ful[(constr.list[[i]]$start): (constr.list[[i]]$end), (constr.list[[i]]$start): (constr.list[[i]]$end)] <- qr.Q(qr(t(constr.list[[i]]$C)), complete=TRUE)			    
	    rm.col <- do.call('c', lapply(constr.list, function(x) x$start))
        ttu <- ttu.ful[, -rm.col]
        modmat.u <- modmat %*% ttu
        penmat.u = crossprod(ttu, penmat %*% ttu)
        nulldim <- nulldim - length(constr.list)
    }
    K = ncol(modmat.u)
    svdP = svd(penmat.u)
    Us = svdP$u[ , 1:(K-nulldim)]
    Un = if (nulldim>0) svdP$u[ , -(1:(K-nulldim))] else rep(0,K)
    X = modmat.u %*% Un
    Z = scale(modmat.u %*% Us, FALSE, sqrt(svdP$d[1:(K-nulldim)]))  # see Wood (2004)
    svdZZ = svd(tcrossprod(Z))
    d = svdZZ$d; d[-(1:ncol(Z))] = 0
  	X. = crossprod(svdZZ$u, X)
    Y. = crossprod(svdZZ$u, Y)
  	rll = Vectorize(function(log.sp) {
            cat(paste("Log smoothing parameter", log.sp),"\n")
    		ev.v = 1 + exp(-log.sp)*d
    		di = diag(1/sqrt(ev.v))
    		X.. = di %*% X.
    		X..X.. = crossprod(X..)
    		Y.. = di %*% Y.
    		m2 = Y.. - X.. %*% solve(X..X.., crossprod(X.., Y..))
    		- (n-ncol(X)) * log(colSums(Y.. * m2)) - sum(log(ev.v)) - log(det(X..X..))
	   })
    tabl = try(rll(lsp)) 
    if (inherits(tabl, "try-error")) {
    	print(svd(modmat.u)$d)
    	stop("If any of the above singular values of the model matrix is (near) zero, error is likely due to rank deficiency of model matrix.")
    }
    colnames(tabl) = lsp
    best.lsp = lsp[as.numeric(apply(tabl, 1, which.max))]
    R = chol(crossprod(modmat.u) + 1e-10 * penmat.u)
    Rinv = solve(R)
    # Rinv = solve(chol(crossprod(modmat.u) + 1e-10 * penmat.u))   # RWC, p. 337
    svdRPR = svd(crossprod(Rinv, penmat.u %*% Rinv))
    RinvU = Rinv %*% svdRPR$u
    tau = svdRPR$d
    if (nulldim > 0) tau[(K-nulldim+1) : K] = 0
    if (nulldim != sum(tau < 1e-10)) warning("Mismatch between null dimension and rank deficiency of transformed penalty matrix")
    A = modmat.u %*% RinvU     # see RWC, p. 336
    M = 1 / (1 + tau %o% exp(best.lsp))
    coef.u = RinvU %*% (M * crossprod(A, Y))
    yhat = modmat.u %*% coef.u
    sigma2 = colSums((Y - yhat)^2)/(n - colSums(M))
    
    # Degrees of freedom per parameter (Wood, 2004)
    # See Jan. 19/12 notes
    edf <- (RinvU * (crossprod(R, svdRPR$u))) %*% M

    fitt = list(fitted = if (store.fitted) yhat else NULL,
                edf = edf, 
                pwdf=colSums(M), 
                pwlsp=best.lsp,
                coef=ttu %*% coef.u,
                reml=if (store.reml) tabl else NULL,  
                modmat = modmat, 
                penmat = penmat,
                # modmat.u = modmat.u, penmat.u = penmat.u, coef.u = coef.u,
                # constr.list = constr.list,
                # Z = Z,
                RinvU = RinvU, tau = tau, 
                sigma2 = sigma2, ttu = ttu)
    class(fitt) = "qplsc.mp"
    fitt
}

