semipar.mp <- function(formula, Y, lsp, data = NULL, range.basis = NULL, knots = "quantile", rm.constr = FALSE, random = NULL, store.reml = FALSE, store.fitted = FALSE) {
    require(stringr)
    tf <- terms.formula(formula, specials="sf")
    trmstrings <- attr(tf, "term.labels")
    # if (!is.null(data)) envir <- environment(formula)
    envir <- environment(formula)
    if (!is.null(range.basis)) envir$range.basis = range.basis
    envir$quantile = quantile
    where.sf <- attr(tf, "specials")$sf
    # where.nsf <- (1:length(trmstrings))[-where.sf]
    where.nsf <- which((1:length(trmstrings))%in%where.sf==FALSE)
    list.all <- vector("list", length(trmstrings))
    
    if (length(where.nsf)!=0) {
        for (i in 1:length(where.nsf)) {
            if (!is.null(data)) {
                list.all[[where.nsf[i]]]$modmat <- as.matrix(eval(as.call(parse(text=trmstrings[where.nsf[i]]))[[1]],envir = data, enclos = envir))} else {
        	    list.all[[where.nsf[i]]]$modmat <- as.matrix(eval(as.call(parse(text=trmstrings[where.nsf[i]]))[[1]], envir = envir))
        	}
            ncol <- ncol(list.all[[where.nsf[i]]]$modmat)
            list.all[[where.nsf[i]]]$penmat <- matrix(0, ncol, ncol) 
            list.all[[where.nsf[i]]]$pen.order <- ncol
        }
    }
    if (length(where.sf)!=0) {
        if (is.null(data)) for (i in 1:length(where.sf)) {
                temptext = trmstrings[where.sf[i]]
                if (!is.null(range.basis)) temptext <- sub("\\)", ", range.basis = range.basis\\)", temptext) 
                list.all[[where.sf[i]]]<-eval(parse(text = temptext), envir = envir)
        	    }
        else for (i in 1:length(where.sf)) {
            formula.term <- sub("sf\\(", "", sub("\\)", "", trmstrings[where.sf[i]]))
            split.term <- as.vector(strsplit(formula.term, ",")[[1]])
            var.x <- split.term[1]
            text <- sub(var.x,paste("data$",var.x,sep=""),trmstrings[where.sf[i]])
            if (!is.na(str_locate(formula.term, ", effect")[1])) { 
    	        effect.ind <- which(is.na(str_locate(split.term, " effect = ")[,1])==FALSE)
                var.effect <- sub(" effect = ", "", split.term[effect.ind])
                text <- sub(var.effect,paste("data$",var.effect,sep=""),text)
            }
            list.all[[where.sf[i]]]<-eval(parse(text = text))
        }
    } 

    modmat <- do.call('cbind', lapply(list.all, function(xx) xx$modmat))
    subpen.dim <- lapply(list.all, function(xx) dim(xx$penmat)[1])
    penmat <- matrix(0, do.call('sum', subpen.dim), do.call('sum', subpen.dim))
    for (i in 1:length(subpen.dim)) {
        start.ind <- cumsum(subpen.dim)[i] - subpen.dim[[i]] + 1
        end.ind <- cumsum(subpen.dim)[i]
        penmat[start.ind:end.ind, start.ind:end.ind] <- list.all[[i]]$penmat
	  list.all[[i]]$start <- start.ind + 1*(attr(tf, "intercept")==1)
	  list.all[[i]]$end <- end.ind + 1*(attr(tf, "intercept")==1)
    }  
    if (!is.null(random)) {
        if (is.matrix(random)) {
            modmat <- cbind(modmat, random)
            penmat <- rbind(cbind(penmat, matrix(0, NCOL(penmat), NCOL(random))), cbind(matrix(0, NCOL(random), NCOL(penmat)), diag(1, NCOL(random))))
        }
        else {
            ran.string <- attr(terms.formula(random), "term.labels")
            ran.terms <- strsplit(ran.string, " \\| ")
            factor <- if (!is.null(data)) as.factor(eval(as.call(parse(text=ran.terms[[1]][2]))[[1]] ,envir = data, enclos = envir))  
                      else as.factor(eval(as.call(parse(text=ran.terms[[1]][2]))[[1]], envir = envir))

          modmat <- cbind(modmat, model.matrix(~ factor-1))
          penmat <- rbind(cbind(penmat, matrix(0, NCOL(penmat), length(levels(factor)))), cbind(matrix(0, length(levels(factor)), NCOL(penmat)), diag(1, length(levels(factor)))))
        }      
    }

    constr.list = NULL
    if (!rm.constr) {
        ind.constr <- do.call('c',lapply(list.all,function(xx) is.null(xx$effect) && !is.null(xx$basis)))
        if (sum(ind.constr)!=0) {
            constr.list <- vector("list", sum(ind.constr))
            for (i in 1:sum(ind.constr)) {
                start <- cumsum(subpen.dim)[which(ind.constr)[i]] - subpen.dim[[which(ind.constr)[i]]] + 1
                end <- cumsum(subpen.dim)[which(ind.constr)[i]]
                constr.list[[i]]$start <- if (attr(tf, "intercept")==1) start+1  else start
                constr.list[[i]]$end <- if (attr(tf, "intercept")==1) end+1  else end
                constr.list[[i]]$C <- matrix(colSums(modmat[,start:end]),1)            
            }
        # if (attr(tf, "intercept")==1) constrmat <- cbind(0, constrmat)
        # else warning("You may wish to adjust constraints for a model without intercept.")
        }
    }
    nulldim <- do.call('sum', lapply(list.all, function(xx) xx$pen.order)) 
    if (attr(tf, "intercept")==1) {
        modmat <- cbind(1, modmat)
        penmat <- rbind(0, cbind(0, penmat))
        nulldim <- nulldim + 1
    }
    qplsc.obj <- qplsc.mp(Y = Y, modmat = modmat, penmat = penmat, constr.list = constr.list, lsp = lsp, nulldim = nulldim, store.reml = store.reml, store.fitted = store.fitted)
  #   if (is.null(random)) {
  #     row.names <- rep(NA, tail(cumsum(subpen.dim),1))
  #     for (i in 1:length(trmstrings)) {
  #         if(subpen.dim[[i]] == 1) row.names[cumsum(subpen.dim)[i]] <- trmstrings[i]
  #         else row.names[(cumsum(subpen.dim)[i]-subpen.dim[[i]]+1):(cumsum(subpen.dim)[i])] <- paste(trmstrings[i], ".", 1:subpen.dim[[i]], sep="")
  #     }
  #     if (attr(tf, "intercept")==1) row.names <- c("(Intercept)", row.names)
  #     rownames(qplsc.obj$coef) <- row.names; rownames(qplsc.obj$edf) <- row.names  
  #   }
    qplsc.obj$where.sf <- where.sf; qplsc.obj$where.nsf <- where.nsf
    qplsc.obj$list.all <- list.all; qplsc.obj$formula <- formula
    qplsc.obj$Y <- Y; qplsc.obj$lsp <- lsp; qplsc.obj$data <- data
    class(qplsc.obj) = c("semipar.mp", "qplsc.mp")
    qplsc.obj
}
