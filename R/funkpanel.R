funkpanel <- function(fkmobj, semiobj, arr4d, predictor, titl="", xlab="", ylab="", ncluster = nrow(fkmobj$centers), 
                      slice=dim(fkmobj$arr.cluster)[3]%/%2, ylim.scatter=NULL, deriv.legend=0, ylim.legend=NULL, 
                      scattermain=NULL, colvec = NULL)  {
    if (!inherits(fkmobj, "funkmeans4d")) stop("'fkmobj' must be an object of class 'funkmeans4d'")
    if (is.null(colvec)) {
    	if (ncluster<=9) colvec = c("dodgerblue", "green", "red", "orange", "yellow", "orchid", "brown", "grey", "purple")[1:ncluster]
        else stop("Please set 'colvec' to a numeric or character vector of colors")
    }
    tkrp = tkrp1 = NULL
    which.smooth = fkmobj$which.smooth
    temp = semiobj$list.all[[which.smooth]]
    screened = min(fkmobj$arr.cluster, na.rm=TRUE) == 0
    if (!is.null(semiobj$incl.inds)) attr(arr4d,"has.data")= semiobj$include
    data.inds = which(attributes(arr4d)$has.data==TRUE, arr.ind=TRUE)
    x.coord = attributes(arr4d)$coord[[1]]
    y.coord = attributes(arr4d)$coord[[2]]
    z.coord = attributes(arr4d)$coord[[3]]
    axis.flag = TRUE
    ttl = "z"
   
    draw <- function(panel) {
        image(x = x.coord, y = y.coord, z = fkmobj$arr.cluster[ , , panel$slice], 
              col = if (screened) c("grey", colvec) else colvec, breaks = .5 + (-screened):ncluster, 
              main = paste(ttl, "=", z.coord[panel$slice]), xlab="x", ylab="y", axes = axis.flag)
        panel
    }
    
    redraw <- function(panel) {
        rp.tkrreplot(panel, tkrp)
        panel
    }
    
    scatterplot <- function(panel,x,y)  {
        arr.ind = c(which.min(abs(x.coord-x)), which.min(abs(y.coord-y)), panel$slice)
        incl.inds = if (is.null(semiobj$incl.inds)) data.inds else semiobj$incl.inds
        which.vox = which(incl.inds[,1]==arr.ind[1] & incl.inds[,2]==arr.ind[2] & incl.inds[,3]==arr.ind[3])
        if (length(which.vox)>0) plot(semiobj, semiobj$Y, which.vox=which.vox, ylab=ylab, main=scattermain, ylim=ylim.scatter)
        else cat("Scatterplots available only for voxels within the displayed clusters!\n")
        return(panel)
   }
 
    legend.draw <- function(panel){
    	fdobj = list()
    	fdobj$coef <- semiobj$coef[(temp$start):(temp$end),]
    	fdobj$basis <- temp$basis
    	plot(fkmobj, fdobj = fdobj, xlab = xlab, ylab = ylab, deriv = deriv.legend, 
    	     ylim = ylim.legend, ncluster = ncluster, colvec = colvec)
    	panel
    }
    
    imgplot <- rp.control(title = titl, slice=slice)
    rp.tkrplot(imgplot, tkrp, draw, action=scatterplot, hscale = 0.8, vscale = 1.2, pos=list(row=0, column=0))
    rp.tkrplot(imgplot, tkrp1, legend.draw, hscale = 0.8, vscale = 1.2, pos=list(row=0, column=1))
    rp.slider(panel = imgplot, var=slice, action = redraw, from = 1, to = dim(fkmobj$arr.cluster)[3], resolution = 1, title=ttl, showvalue=is.null(attributes(arr4d)$coord), pos=list(row=1, column=1))
}

