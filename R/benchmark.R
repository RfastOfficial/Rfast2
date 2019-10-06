

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
    colnames(res)<-c("min","mean","max")
    class(res)<-"benchmark"
    res
}

print.benchmark<-function(x,...){
    class(x)<-NULL
    min_max_res<-x[,c("min","max")]
    unit<-NULL
    unit_number<-NULL
    if(any(min_max_res>=1)){
        unit<-"seconds"
        unit_number<-1 #nothing
    }else if(any(min_max_res>=10^-2) && any(min_max_res<1)){
        unit<-"milliseconds"
        unit_number<-1000 # 10^3
    }else if(any(min_max_res>=10^-3) && any(min_max_res<10^-2)){
        unit<-"milliseconds"
        unit_number<-10000 # 10^4
    }else if(any(min_max_res>=10^-6) && any(min_max_res<10^-3)){
        unit<-"microseconds"
        unit_number<-1000000 #10^6
    }else{
        stop("Error: times is weird...\n")
    }
    
    x<-x*unit_number
    cat("  ",unit,"\n")
    print(x)
}
