# Y: n x V outcome matrix (V=number of voxels, connections, etc.)
# formula: object of form "~ x1 + x2" (quotation marks not needed)
lm.mp <- function(Y,formula, store.fitted=FALSE) {
    ## Y = eval(formula[[2]], parent.frame())
    X = model.matrix(formula)
    n = dim(X)[1]
    p = dim(X)[2]
    XtX.inv = solve(crossprod(X))
    I.H = diag(n) - X %*% tcrossprod(XtX.inv, X)
    coef = XtX.inv %*% crossprod(X, Y)
    sigma2 = apply(I.H %*% Y, 2, crossprod) / (n-p)
    se.coef = sqrt(diag(XtX.inv) %o% sigma2)
    fitted=if (store.fitted) X %*% coef else NULL
    otpt = list(coef=coef, sigma2=sigma2, se.coef=se.coef, X=X, fitted=fitted) 
    class(otpt) = "lm.mp"
    otpt
}

