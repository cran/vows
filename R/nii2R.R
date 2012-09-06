nii2R <-
function(niifilename, which.vols=NULL, savename=NULL, remove.zero=TRUE, maskname=NULL, ind=NULL, ind.auto=TRUE, coord=NULL)   {
    require(Rniftilib)
    bigobj = nifti.image.read(niifilename)
	if (length(dim(bigobj))==3) {
		tmp = array(dim=c(dim(bigobj), 1))
		tmp[ , , , 1] = bigobj[]
		bigobj = tmp
		rm(tmp)
	}
    
    if (is.null(ind) & ind.auto)    ind = get.ind(bigobj[,,,1])
    else if (is.null(ind) & !ind.auto)    for (i in 1:3)  ind[[i]]=1:dim(bigobj[])[i]
    x.ind = ind[[1]]
    y.ind = ind[[2]]
    z.ind = ind[[3]]
    if (is.null(which.vols)) which.vols = 1:(dim(bigobj)[4])
    
    if (is.null(coord)) {
  		coord = list()
  		for (l in 1:3) coord[[l]] = 1:(dim(bigobj[])[l])
    }
    
    coord[[1]] = coord[[1]][x.ind]
    coord[[2]] = coord[[2]][y.ind]
    coord[[3]] = coord[[3]][z.ind]
        
    dir.create("./temp")
    cat("Saving each image to temp directory...\n")
    for (l in which.vols)   {
        cat("Image #",l,"\n")
        temp = bigobj[x.ind, y.ind, z.ind, l]
        save(temp, file=paste("./temp/Sep-", l, ".RData", sep=""))
    }
    dim.nii = dim(bigobj[,,,1])
    rm(temp, bigobj)
    d4 = array(NA, c(length(x.ind), length(y.ind), length(z.ind), length(which.vols)))
    cat("Adding each image 4D array...\n")
    i=0
    for (l in which.vols)   {
        i = i+1
        cat("Image #",i,"\n")
        load(paste("./temp/Sep-", l, ".RData", sep=""))
        d4[,,,i]=temp
    }
    
    if (!is.null(maskname))   has.data = nifti.image.read(maskname)[x.ind,y.ind,z.ind]!=0
    else if (remove.zero)   has.data =  apply(d4, 1:3, function(mat) !all(mat==0|is.infinite(mat)|is.na(mat)))
         else   has.data = array(TRUE, dim(d4)[1:3])
    
    attr(d4, "x.ind") = x.ind
    attr(d4, "y.ind") = y.ind
    attr(d4, "z.ind") = z.ind
    attr(d4, "which.vols") = which.vols
    attr(d4, "dim.nii") = dim.nii
    attr(d4, "coord") = coord
    attr(d4, "has.data") = has.data
    if (!is.null(savename)) save(d4, file=paste(savename,".RData",sep=""))
    return(d4)
}

