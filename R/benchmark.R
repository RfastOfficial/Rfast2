#[export]
benchmark<-function(...,times,envir=parent.frame(),order=NULL){
    exprs<-as.list(match.call(expand.dots = FALSE)$...)
    nm<-names(exprs)
    index_to_empty_names<-which(nm=="")
    if(is.null(nm))# if none of the expression has name
        nm<-sapply(exprs,function(x){paste(deparse(x),collapse = "")})
    else{
        if(length(index_to_empty_names)>0)# get only the expressions that do not have a name
            nm[index_to_empty_names]<-sapply(exprs[index_to_empty_names],function(x){paste(deparse(x),collapse = "")})
    }
    if(is.null(order))
        order<-sample(length(exprs),length(exprs))
    else if(length(order)!=length(exprs))
        stop("Error: order length must be ",length(exprs))
    res<-.Call(Rfast2_benchmark,exprs,envir,times,order)
    rownames(res)<-nm
    colnames(res)<-c("min","mean","max","CPU","memory")
    class(res)<-"benchmark"
    res
}


#[export s3]
print.benchmark<-function(x,...){
    class(x)<-NULL
    x <- as.data.frame(x)
    min_mean_max<-c("min","mean","max")
    min_mean_max_res<-x[,min_mean_max]
    unit<-NULL
    unit_number<-NULL
    if(any(min_mean_max_res>=1)){
        unit<-"seconds"
        unit_number<-1 #nothing
    }else if(any(min_mean_max_res>=10^-2) && any(min_mean_max_res<1)){
        unit<-"milliseconds"
        unit_number<-1000 # 10^3
    }else if(any(min_mean_max_res>=10^-3) && any(min_mean_max_res<10^-2)){
        unit<-"milliseconds"
        unit_number<-10000 # 10^4
    }else if(any(min_mean_max_res>=10^-6) && any(min_mean_max_res<10^-3)){
        unit<-"microseconds"
        unit_number<-1000000 #10^6
    }else{
        stop("Error: times is weird...\n")
    }

    x[,min_mean_max]<-x[,min_mean_max]*unit_number
    x[,"CPU"] <- paste(round(x[,"CPU"],2),"%",sep = "")

    memories <- c("B"=1,"KB"=1024,"MB"=1024000,"GB"=1024000000)
    memory <- x[,"memory"] / memories
    memory_index <- which(memory >= 0 & memory <= 100)[1] #if zero the result will be more than one so pick up the first (Byte)
    memory <- memory[memory_index]
    x[,"memory"] <- paste(round(memory,2),names(memory),sep = "")
    


    cat("  ",unit,"\n")
    print(x)
}
