#[export]
colTrimMean<-function(x,a=0.05,parallel=FALSE){
	.Call(Rfast2_colTrimMean,x,a,parallel)
}
