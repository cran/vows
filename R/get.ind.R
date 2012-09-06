get.ind <-
function(arr3d)   {
    flag = rep(NA, dim(arr3d)[1])
    for (i in  1:length(flag))    {
        cat("First-dimension slice",i,"\n")
        flag[i] = any(arr3d[i,,]!=0)
    }
    x.ind = which(flag)
    
    flag = rep(NA, dim(arr3d)[2])
    for (i in  1:length(flag))    {
        cat("Second-dimension slice",i,"\n")
        flag[i] = any(arr3d[,i,]!=0)
    }
    y.ind = which(flag)
	
    flag = rep(NA, dim(arr3d)[3])
    for (i in  1:length(flag))    {
        cat("Third-dimension slice",i,"\n")
        flag[i] = any(arr3d[,,i]!=0)
    }
    z.ind = which(flag)
    return(list(x.ind, y.ind, z.ind))
}

