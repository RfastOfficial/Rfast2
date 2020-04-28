#[export]
rowTrimMean<-function(x,a=0.05,parallel=FALSE){
	.Call(Rfast2_rowTrimMean,x,a,parallel)
}
