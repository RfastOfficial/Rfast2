#[export]
colQuantile<-function(x,probs,parallel=FALSE){
	.Call(Rfast2_col_Quantile,x,probs,parallel)
}
