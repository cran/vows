# TODO: p-values for factors, etc.
summary.lm.mp <- function(object, ...) {
	X = object$X
	coef = object$coef
	se.coef = object$se.coef
    n = dim(X)[1]
    p = dim(X)[2]
    tstat = coef[-1,] / se.coef[-1,]  
    pvalue = 2 * pt(-abs(tstat), n-p)  
	sigma2.mle = object$sigma2 * (n-p) / n
    aicc = n*log(sigma2.mle) + n*(n+p)/(n-p-2)	    
	list(tstat=tstat, pvalue=pvalue, aicc=aicc)
}

