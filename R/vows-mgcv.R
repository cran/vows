# Adapted from smooth.construct.ps.smooth.spec
# B-splines with knots at quantiles of data
smooth.construct.bq.smooth.spec <- function(object, data, knots) {
    require(fda)
    if (length(object$p.order) == 1) 
        m <- rep(object$p.order, 2)
    else m <- object$p.order
    m[is.na(m)] <- 2
    object$p.order <- m
    if (object$bs.dim < 0) 
        object$bs.dim <- max(10, m[1] + 1)
    if (object$bs.dim <= m[1]) 
        stop("basis dimension too small for B-spline order")
    if (length(object$term) != 1) 
        stop("Basis only handles 1D smooths")
    x <- data[[object$term]]
    k <- knots[[object$term]]
    if (is.null(k)) {
	     bsb <- create.bspline.basis(range(x), norder=m[1]+2, nbasis = object$bs.dim, breaks = quantile(x, seq(0,1,length.out=object$bs.dim-m[1])))
	     k <- c(bsb$rangeval[1], bsb$params, bsb$rangeval[2])
    }            
	else bsb <- create.bspline.basis(range(x), norder=m[1]+2, breaks=k)

    object$X <- eval.basis(x, bsb)
    object$S <- list(getbasispenalty(bsb, m[2]))
    object$rank <- object$bs.dim - m[2]
    object$null.space.dim <- m[2]
    object$knots <- k
    object$m <- m
    class(object) <- "bspline.smooth"
    object
}

# B-splines with equally-spaced knots
smooth.construct.be.smooth.spec <- function(object, data, knots) {
    require(fda)
    if (length(object$p.order) == 1) 
        m <- rep(object$p.order, 2)
    else m <- object$p.order
    m[is.na(m)] <- 2
    object$p.order <- m
    if (object$bs.dim < 0) 
        object$bs.dim <- max(10, m[1] + 1)
    if (object$bs.dim <= m[1]) 
        stop("basis dimension too small for B-spline order")
    if (length(object$term) != 1) 
        stop("Basis only handles 1D smooths")
    x <- data[[object$term]]
    k <- knots[[object$term]]
    if (is.null(k)) {
	     bsb <- create.bspline.basis(range(x), norder=m[1]+2, nbasis = object$bs.dim)
	     k <- c(bsb$rangeval[1], bsb$params, bsb$rangeval[2])
    }            
	else bsb <- create.bspline.basis(range(x), norder=m[1]+2, breaks=k)

    object$X <- eval.basis(x, bsb)
    object$S <- list(getbasispenalty(bsb, m[2]))
    object$rank <- object$bs.dim - m[2]
    object$null.space.dim <- m[2]
    object$knots <- k
    object$m <- m
    class(object) <- "bspline.smooth"
    object
}

Predict.matrix.bspline.smooth<-function (object, data) {
    require(fda)
    x = data[[object$term]]
    k = object$knots
    m = object$m
    bsb <- create.bspline.basis(range(x), norder=m[1]+2, breaks=k)
    X <- eval.basis(x, bsb)
    X
}