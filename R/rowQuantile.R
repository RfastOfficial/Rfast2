#[export]
rowQuantile<-function(x,probs,parallel=FALSE){
	.Call(Rfast2_row_Quantile,x,probs,parallel)
}
