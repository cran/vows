plot.rlrt4d <-
function(x, array4d, disp = c("stat", "p", "fdr", "pwdf"), titl=NULL, slices = NULL, colbar = TRUE, col.image = femmecol(100)[100:1], neglog10=FALSE, threshold=NULL, mar=c(2,2,2,2), digit=2, nrow=NULL, ...)  {
    disp = match.arg(disp)
    x.ind = attributes(array4d)$x.ind
    y.ind = attributes(array4d)$y.ind
    z.ind = attributes(array4d)$z.ind
    coord = attributes(array4d)$coord
    has.data = attributes(array4d)$has.data
    
    x.coord = coord[[1]]
    y.coord = coord[[2]]
    z.coord = coord[[3]]
    axis.flag = TRUE
    ttl = "z ="
    xlb="x"; ylb="y"
    
    if (!(disp %in% c("stat", "p", "fdr"))) stop("You must choose RLRT statistics, p-value or fdr to display!")
    arr = array(NA, dim=dim(has.data))
    if (disp=="stat") arr[has.data] = x$stat
    if (disp=="p") arr[has.data] = x$p
    if (disp=="fdr") arr[has.data] = x$fdr
    if (neglog10) {
    	arr[has.data] = -log10(arr[has.data])
        if (disp=="stat") warning("Do you really want to take the negative base-10 log of the RLR statistic?")
    }
    
    arr.le.th = arr.ge.th = arr
    if (!is.null(threshold)) {
        arr.le.th[!is.na(arr) & (arr>threshold)] = NA
        arr.ge.th[!is.na(arr) & (arr<=threshold)] = NA 
        arr.ge.th[!is.na(arr) & (arr>threshold)] = 1
    }
    zlim = range(c(range(arr.le.th, na.rm=TRUE), threshold), na.rm = TRUE)  
    
    if (is.null(slices)) slices = round(seq(5,dim(arr)[3]-4,,11))
    if (is.null(nrow)) nrow = ceiling(sqrt(length(slices)+colbar))
    ncol = ceiling((length(slices)+colbar)/nrow)
    par(mfrow=c(nrow, ncol), mar = mar)
    
    for (i in 1:length(slices)){
        image(x = x.coord, y=y.coord, z=arr.le.th[ , , slices[i]],
              col=col.image, main= ifelse(is.null(titl), "", paste(ttl, z.coord[slices[i]])), 
              xlab=xlb, ylab=ylb, zlim=zlim, axes=axis.flag, ...)
        image(x = x.coord, y=y.coord, z=arr.ge.th[ , , slices[i]], col="grey", add=TRUE)    
    }
    if (colbar) {
        emptyplot(main="    ")
        colorlegend(posx=c(0.6,0.7), col=col.image,
                    zlim=zlim, zval = seq(min(zlim), max(zlim), length.out=5),main="", left=FALSE, digit=digit)
    }
}

