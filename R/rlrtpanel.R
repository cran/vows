rlrtpanel <-
function(rlrtobj, array4d, predictor, disp = c("stat","p","fdr"), titl="", xlab="Age", ylab="", slice=(dim(array4d)[3]%/%2), ylim.scatter=NULL, col.image = femmecol(100), neglog10=FALSE, threshold=NULL)  {
    require(rpanel)
    require(shape)
    tkrp = tkrp1 = NULL
    disp = match.arg(disp)
    coord = attributes(array4d)$coord
    has.data = attributes(array4d)$has.data
    
    x.coord = coord[[1]]
    y.coord = coord[[2]]
    z.coord = coord[[3]]
    axis.flag = TRUE
    ttl = "z ="
    xlb="x"; ylb="y"
    if (!(disp %in% c("stat", "p", "fdr"))) stop("You must choose rlrt statistics, p-value or fdr to display!")
    arr = array(NA, dim=dim(has.data))
    if (disp=="stat") arr[has.data] = rlrtobj$stat
    if (disp=="p") arr[has.data] = rlrtobj$p
    if (disp=="fdr") arr[has.data] = rlrtobj$fdr
    if (neglog10) {
    	arr[has.data] = -log10(arr[has.data])
        if (disp=="stat") warning("Do you really want to take the negative base-10 log of the RLR statistic?")
    }
    if (!is.null(threshold))    arr[!is.na(arr) & (arr>threshold)] = threshold
    zlim = range(arr[arr<Inf], na.rm=TRUE)    
    
    draw <- function(panel) {
        with(panel,image(x = x.coord, y=y.coord, z=arr[ , , slice],
              col=col.image, main=paste(ttl, z.coord[slice]), xlab=xlb, ylab=ylb, zlim=zlim, axes=axis.flag))
        panel
    }
    
    redraw <- function(panel) {   
        rp.tkrreplot(panel, tkrp)
        panel
    }
    
    scatterplot <- function(panel,x,y)  {
      with(panel, plot(predictor, array4d[which.min(abs(x.coord-x)), which.min(abs(y.coord-y)), slice, ],
           xlab=xlab, ylab=ylab, ylim=ylim.scatter, main=paste("Coordinates", round(x), round(y), z.coord[slice])))
      return(panel)
   }
   
    legend.draw <- function(panel) {
         emptyplot(main="    ")
         colorlegend(posx=c(0.6,0.7), col=col.image,
               zlim=zlim, zval = seq(min(zlim), max(zlim), length.out=5),main="", left=FALSE, digit=2)
         panel
    }
    
    imgplot <- rp.control(title = titl, slice=slice)
    rp.tkrplot(imgplot, tkrp, draw, action=scatterplot, hscale = 0.8, vscale = 1.2, pos=list(row=0, column=0))
    rp.tkrplot(imgplot, tkrp1, legend.draw, hscale = 0.8, vscale = 1.2, pos=list(row=0, column=1))
    rp.slider(panel = imgplot, variable=slice, action = redraw, from = 1, to = dim(arr)[3], resolution = 1, title=ttl, showvalue=is.null(coord), pos=list(row=1, column=1))
}

