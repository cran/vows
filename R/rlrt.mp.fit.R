rlrt.mp.fit <-
function(Y, X, Z, loginvsp, evalarg=NULL, get.df=FALSE) {
    require(RLRsim)
    sim = RLRTSim(X, Z, qrX=qr(X), sqrt.Sigma=diag(NCOL(Z)))
    try(require(fda))
	if (!exists("is.fd")) is.fd = function(mm) FALSE
	n = NROW(X); p = NCOL(X); q = NCOL(Z)
	svdZZ = svd(tcrossprod(Z))
	d = svdZZ$d; d[-(1:q)] = 0
	X. = crossprod(svdZZ$u, X)
	X.X. = crossprod(X.); detX.X. = det(X.X.)
	Y. = crossprod(svdZZ$u, if (is.fd(Y)) t(Y$coef) else Y)
	if (is.fd(Y)) {
		if (is.null(evalarg)) stop("'evalarg' must be specified when 'Y' is of class 'fd'")
		B = eval.basis(evalarg, Y$basis) 
	}
	m1 = Y. - X. %*% solve(X.X., crossprod(X., Y.))
	if (is.fd(Y)) term1 = (n-p) * log(colSums(tcrossprod(Y.,B) * tcrossprod(m1,B)))
	else term1 = (n-p) * log(colSums(Y. * m1))
	rldiff2 = Vectorize(function(log.inv.sp) {
		ev.v = 1 + exp(log.inv.sp)*d
		di = diag(1/sqrt(ev.v))
		X.. = di %*% X.
		X..X.. = crossprod(X..)
		Y.. = di %*% Y.
		m2 = Y.. - X.. %*% solve(X..X.., crossprod(X.., Y..))
		if (is.fd(Y)) term2 = (n-p) * log(colSums(tcrossprod(Y..,B) * tcrossprod(m2,B)))
	    else term2 = (n-p) * log(colSums(Y.. * m2))
		term1 - term2 - sum(log(ev.v)) - log(det(X..X..)) + log(detX.X.)
	})
	if (get.df) {  # see RWC, p. 336
		Rinv = solve(chol(crossprod(cbind(X, Z))))
		Rinv22 = Rinv[-(1:p), -(1:p)]
		svec = c(svd(Rinv22)$d^2, rep(0, p))
	}
	tabl = matrix(rldiff2(loginvsp), ncol = length(loginvsp))
	colnames(tabl) = loginvsp
	logsp = -loginvsp[apply(tabl, 1, which.max)]
    stat = apply(tabl, 1, max)
    pvalue = stat
    for (i in 1:length(stat)) pvalue[i] = mean(stat[i] < sim)
    fdr = p.adjust(pvalue, method="BH")
	list(table = tabl, stat = stat, logsp = logsp, df = if (get.df) rowSums(1/(1+exp(logsp) %o% svec)) else NULL, sim = sim, pvalue=pvalue, fdr=fdr)
}

