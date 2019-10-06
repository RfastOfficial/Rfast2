

#[export]
is.skew.symmetric<-function(x){
	.Call(Rfast2_is_skew_symmetric,x)
}
