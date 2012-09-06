permF.mp <-
function(formula, nperm=499, alpha=.05, report.every=50)  {
    Y = eval(formula[[2]], parent.frame())
    X = model.matrix(formula)
    n = dim(X)[1]; p = dim(X)[2]
    I.H = diag(n) - X %*% solve(crossprod(X), t(X))
    rss0 = apply(scale(Y, TRUE, FALSE), 2, crossprod)
    maxF.perm = c()
    for (i in 1:nperm)  {
        if (i %% report.every==0) cat("Permutation", i, "\n")
        Y.perm = Y[sample(1:n), ]
        rss1.perm = apply(I.H %*% Y.perm, 2, crossprod)
        F.perm = ((n-p)/(p-1)) * (rss0-rss1.perm) / rss1.perm
        maxF.perm[i] = max(F.perm)
    }
    rss1 = apply(I.H %*% Y, 2, crossprod)
    F.obs = ((n-p)/(p-1)) * (rss0-rss1) / rss1
    threshold = quantile(maxF.perm, 1-alpha)
    pvalue = (rowSums(outer(F.obs, maxF.perm, function(x,y) x<y)) + 1) / (nperm + 1)
    return(list(maxF.perm=maxF.perm, F.obs=F.obs, threshold=threshold, pvalue=pvalue))
}

