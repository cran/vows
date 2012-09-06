F.mp <-
function(formula, which)    {
    Y = eval(formula[[2]], parent.frame())
    X = model.matrix(formula)
    n = dim(X)[1]
    p = dim(X)[2]
    df1 = length(which)
    XtX.inv = solve(crossprod(X))
    I.H = diag(n) - X %*% tcrossprod(XtX.inv, X)
    rss1 = apply(I.H %*% Y, 2, crossprod)
    X0 = matrix(X[ , -which], nrow = nrow(X))
    XtX.inv0 = solve(crossprod(X0))
    I.H0 = diag(n) - X0 %*% tcrossprod(XtX.inv0, X0)
    rss0 = apply(I.H0 %*% Y, 2, crossprod)
    F = ((n-p)/df1) * (rss0-rss1) / rss1
    otpt = list(F=F, df1=df1, df2=n-p, pvalue=pf(F, df1, n-p, lower.tail=FALSE), X=X)
    otpt
}

