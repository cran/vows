lm4d <- function(arr4d, formula, store.fitted=FALSE) {
    dim.4d = dim(arr4d)
    has.data<-attributes(arr4d)$has.data
    N=NROW(model.matrix(formula))
    dim(arr4d)=c(prod(dim(arr4d)[1:3]),N)
    Y.d = t(arr4d[(which(has.data)),])  
    rm(arr4d)  # i.e., remove from current environment
    Y.fit=lm.mp(Y=Y.d,formula=formula, store.fitted=store.fitted)
    if (store.fitted) {
        fit.value=array(NA,dim=dim.4d)
        for (j in 1:N) fit.value[,,,j][which(has.data, arr.ind=TRUE)]=Y.fit$fitted[j,]
        Y.fit$fitted=fit.value
    }
    Y.fit$call = match.call()
    return(Y.fit)
}
