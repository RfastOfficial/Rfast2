#[export]
colGroup<-function(x,ina,method="sum",names=TRUE, std = FALSE){
	unique_ina<-unique(ina)
	
	if(method=="var"){
	    m <- rowsum(x, ina)
		m2 <- rowsum(x^2, ina)
		ni <- tabulate(ina)
		ni <- ni[ni > 0]
		y <- (m2 - m^2/ni) / (ni - 1)
		if ( std ) y <- sqrt(y)
	}
	else y <- .Call(Rfast2_col_group,x,ina,length(unique_ina),method)
	
	if(names){
		rownames(y)<-as.character(Rfast::Sort(unique_ina))
	}
	
	y
}
