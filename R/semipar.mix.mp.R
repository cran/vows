semipar.mix.mp <- function(Y, x, param=NULL, random, data.ran, k = 10, norder = 4, pen.order=2, knots = "quantile", store.gamm4 = FALSE) {
    require(gamm4)	 
    if (store.gamm4) {
    	gamm4.list <- vector("list", NCOL(Y)) 
    }
    if (knots == "quantile") bs = "bq"
    else if (knots == "equispaced") bs = "be"
    coefmat <- matrix(NA, k + NCOL(param)-is.null(param), NCOL(Y))
    bsplinecoefmat <- matrix(NA, k, NCOL(Y))
    pwdf <- pwlsp <- rep(NA, NCOL(Y))
    ss = smoothCon(s(x, k = k, bs=bs, m=c(norder-2,pen.order)), data.frame(x = x), knots=NULL)
    B <- ss[[1]]$X 
    C <- ss[[1]]$C
    # all.equal(C, matrix(colSums(ss[[1]]$X), 1))  # TRUE
    Z <- qr.Q(qr(t(C)), complete=TRUE)[ , -1]   
    # Note that Z, as defined above, transforms the non-intercept coeffs 
    # to the original parametrization.
    basis <-  create.bspline.basis(range(x), norder = norder, nbasis = k)

    if(!is.null(param)) {
	   frame <- vector("list", length(all.vars(random)) + 3)
	   names(frame) <- c("y", "x", "param", all.vars(random))
         frame$x <- x
	   frame$param <- param
	   for(i in 1:length(all.vars(random))) frame[[3+i]] <- data.ran[, all.vars(random)[i]]
    }
    if (is.null(param)) {
         frame <- vector("list", length(all.vars(random)) + 2)
	   names(frame) <- c("y", "x", all.vars(random))
         frame$x <- x
	   for(i in 1:length(all.vars(random))) frame[[2+i]] <- data.ran[, all.vars(random)[i]]
    }

    for (j in 1:NCOL(Y)) {
    	if (!j%%40) cat("Vertex", j, "\n")
        if (is.null(param)) {
        	  frame$y <- Y[, j]
 	        gamm4.obj <-gamm4(y ~ s(x, k = k, bs=bs, m=c(norder-2,pen.order)), random=random, data=frame)
 	    }
 	    else {
   	       frame$y <- Y[, j]
             gamm4.obj <-gamm4(y ~ param + s(x, k = k, bs=bs, m=c(norder-2,pen.order)), random=random, data = frame)   
	 }
        coefmat[, j] <- gamm4.obj$gam$coefficients
        bsplinecoefmat[, j] <- coefmat[1, j] * rep(1, k) + Z %*% if (is.null(param)) coefmat[-1, j] else coefmat[-(1:(1+NCOL(param))), j]
        pwdf[j] <- sum(gamm4.obj$gam$edf)
        pwlsp[j] <- log(gamm4.obj$gam$sp)
        if (store.gamm4) gamm4.list[[j]] <- gamm4.obj
    } 
    list(coef = coefmat, bsplinecoef = bsplinecoefmat, pwdf = pwdf, pwlsp = pwlsp, B = B, C = C, 
         Z = Z, basis = basis, gamm4.list = if(store.gamm4) gamm4.list else NULL)
}
