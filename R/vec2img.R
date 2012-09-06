vec2img<-function(vec,logicarray){
    if(length(vec)!=sum(logicarray)) stop("dimension not fit")
    trueinds<-which(logicarray==TRUE,arr.ind=TRUE)
    result<-array(NA,dim(logicarray))
    result[logicarray]<-vec
    return(result)
}
