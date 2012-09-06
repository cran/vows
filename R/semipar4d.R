semipar4d <-
function(arr4d, formula, lsp, data, range.basis = NULL, knots = "quantile", rm.constr = FALSE, random = NULL, store.reml = FALSE, store.fitted = FALSE)   {
    dim.4d = dim(arr4d)
    has.data <- attributes(arr4d)$has.data
    N = NROW(data)
    dim(arr4d) = c(prod(dim(arr4d)[1:3]), N)
    Y.d = t(arr4d[(as.vector(has.data)), ])
    rm(arr4d) # i.e., remove from current environment  
    
    Y.fit = semipar.mp(formula = formula, Y = Y.d, lsp = lsp, data = data, range.basis = range.basis, knots = knots, rm.constr = rm.constr, store.reml = store.reml, store.fitted = store.fitted)
    if (store.fitted) {
        fit.value = array(NA, dim=dim.4d)
        for (j in 1:N) fit.value[ , , ,j][has.data] = Y.fit$fitted[j,]
        Y.fit$fitted = fit.value
    }
    Y.fit$call = match.call()
    return(Y.fit)
}
 