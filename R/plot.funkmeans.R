plot.funkmeans = function(x, fdobj, deriv=0, ylim=NULL, ncluster = nrow(x$centers), mfrow = NULL,
            colvec = NULL,...) {
    if (!is.null(mfrow)) par(mfrow=mfrow)
    else {
    	nro = ceiling(sqrt(ncluster))
    	nco = ceiling(ncluster/nro)
    	par(mfrow=c(nro,nco))
    }
    if (is.null(colvec)) {
    	if (ncluster<=9) colvec = c("dodgerblue", "green", "red", "orange", "yellow", "orchid", "brown", "grey", "purple")[1:ncluster]
        else stop("Please set 'colvec' to a numeric or character vector of colors")
    }
    for (ii in 1:ncluster) {
        ee = sample(which(x$cluster==ii), min(30, sum(x$cluster==ii)))
        ee.fd = deriv.fd(fd(coef = fdobj$coef[ , ee], basisobj = fdobj$basis), deriv)
        plot(ee.fd, lty=1, col=colvec[x$cluster[ee]], ylim=ylim, main=sum(x$cluster==ii), ...)           
        idx = which(x$cluster==ii)
        mean.cluster = deriv.fd(fd(coef = apply(fdobj$coef[ , idx],1,mean), basisobj = fdobj$basis), deriv)
        lines(mean.cluster, lty=1, col="black", lwd=2)
    }
}
