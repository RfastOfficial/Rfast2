
#[export]
Quantile<-function(x,probs){
	.Call(Rfast2_Quantile,x,probs)
}
