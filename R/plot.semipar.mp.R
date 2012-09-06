plot.semipar.mp <-		
function(x, Y, arr.ind = NULL, which.vox = NULL, which.smooth = NULL, coverage = 0.95, length.new = 100, ...) {		
	if (is.null(arr.ind) & is.null(which.vox)) stop("Must specify either 'arr.ind' or 'which.vox'")	
		
    tf <- terms.formula(x$formula, specials = "sf")		
    trmstrings <- attr(tf, "term.labels")		
    terms <- rep(NA, length(x$where.sf)); smooth.terms <- terms		
    for (i in 1:length(x$where.sf)) {		
	    formula.term <- sub("sf\\(", "", sub("\\)", "", trmstrings[x$where.sf[i]]))	
        split.term <- as.vector(strsplit(formula.term, ",")[[1]])		
        terms[i] <- split.term[1]		
        smooth.terms[i] <- paste("sf(", split.term[1], ")", sep="")		
	}	
		
    alpha = 1-coverage; mul = qnorm(1-alpha/2)		
    intercept <- attr(terms.formula(x$formula, specials="sf"), "intercept")		
		
    ################################################################################		
 		
    B.list <- plot.list <- vector("list", length(x$where.sf))		
    for  (i in 1:length(x$where.sf)) {		
	 xarg <- x$list.all[[x$where.sf[i]]]$argvals	
       x.new <- seq(range(xarg)[1], range(xarg)[2], length.out = length.new)		
	 plot.list[[i]]$x <- xarg	
    	 plot.list[[i]]$x.new <- x.new	
       B.list[[i]] = eval.basis(plot.list[[i]]$x.new, x$list.all[[x$where.sf[i]]]$basis)           		
    }		
     pred.ind <- rep(NA, length(x$where.sf)+length(x$where.nsf))		
     for(i in (x$where.sf)) pred.ind[i] <- x$list.all[[i]]$k		
     for(i in (x$where.nsf)) pred.ind[i] <- ncol(x$list.all[[i]]$modmat)		
     start.ind <- c(1, cumsum(pred.ind)[-length(pred.ind)]+1) + intercept		
     end.ind <- start.ind + pred.ind - 1   		
     modmat.tmp <- matrix(NA, length.new, NCOL(x$modmat))		
     if(intercept) modmat.tmp[, 1] <- 1		
     if(length(x$where.nsf)>0) {		
     for(i in 1:length(x$where.nsf)) {		
    		modmat.tmp[, (start.ind[x$where.nsf[i]]):(end.ind[x$where.nsf[i]])] <- matrix(apply(as.matrix(x$modmat[, (start.ind[x$where.nsf[i]]):(end.ind[x$where.nsf[i]])]), 2, mean), length.new, byrow = TRUE)
     }}		
     for(i in 1:length(x$where.sf)) {		
    		modmat.tmp[, (start.ind[x$where.sf[i]]):(end.ind[x$where.sf[i]])] <- 0
     }		
     Xarray <- array(modmat.tmp, dim = c(dim(modmat.tmp), length(x$where.sf))) 		
     for(i in 1:length(x$where.sf)) {		
        Xarray[,(start.ind[x$where.sf[i]]):(end.ind[x$where.sf[i]]),i] <- B.list[[i]]		
    } 		
    ################################################################################		
       coef <- x$coef[, which.vox]		
          for  (i in 1:length(x$where.sf)) {		
		coef.seq <- (start.ind[x$where.sf[i]]):(end.ind[x$where.sf[i]])
            plot.list[[i]]$y.hat <- Xarray[,coef.seq,i]%*%coef[coef.seq]		
          }		
     if (is.null(which.smooth)) plot.ind <- 1:length(x$where.sf)		
	else plot.ind <- which.smooth	
        par(mfrow=c(1, length(plot.ind)))		
        for (i in plot.ind) {    		
		B <- Xarray[,,i]%*%x$ttu      
		pwvar <- x$sigma2[which.vox] * (B %*% x$RinvU)^2 %*% (1 / (1 + exp(x$pwlsp[which.vox]) * x$tau))
		
		ord <- order(plot.list[[i]]$x.new)
            lower <- min(plot.list[[i]]$y.hat-mul*sqrt(pwvar))		
            upper <- max(plot.list[[i]]$y.hat+mul*sqrt(pwvar))		
            range <- upper-lower		
            plot(plot.list[[i]]$x.new[ord], plot.list[[i]]$y.hat[ord], xlab = terms[i], ylab = smooth.terms[i], type="l", ylim=c(lower-range*0.1, upper+range*0.1), ...)		
            rug(plot.list[[i]]$x)		
            lines(plot.list[[i]]$x.new[ord], plot.list[[i]]$y.hat[ord]-mul*sqrt(pwvar[ord]), lty=2)		
            lines(plot.list[[i]]$x.new[ord], plot.list[[i]]$y.hat[ord]+mul*sqrt(pwvar[ord]), lty=2)    		
        } 		
}		
		
