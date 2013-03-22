funkmeans = function(fdobj, deriv = 1, lambda = 0, ncomp, centers, nstart = 10) {   
    deriv.fdobj = deriv.fd(fdobj, deriv)
    harmfdPar = fdPar(deriv.fdobj)
    harmfdPar$lambda = lambda
    fpca.obj = pca.fd(deriv.fdobj, nharm = ncomp, harmfdPar)
    km.obj = kmeans(fpca.obj$scores, centers = centers, nstart = nstart, iter.max=100)
    km.obj$basis = fdobj$basis
    km.obj$coef = fdobj$coef
    km.obj$fpca = fpca.obj
    km.obj$R2 = (1 - km.obj$tot.withinss / km.obj$totss) * sum(km.obj$fpca$varprop) 
    class(km.obj) = "funkmeans"
    km.obj
}
