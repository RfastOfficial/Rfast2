#[export]
trim.mean<-function(x,a=0.05){
	.Call(Rfast2_trimmean,x,a)
}
