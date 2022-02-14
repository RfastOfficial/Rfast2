#[export]
is.lower.tri <- function(x,diag = FALSE) {
	.Call(Rfast2_is_lower_tri, x, diag)
}



#[export]
is.upper.tri <- function(x,diag = FALSE) {
	.Call(Rfast2_is_upper_tri, x, diag)
}



#[export]
is.skew.symmetric<-function(x){
	.Call(Rfast2_is_skew_symmetric, x)
}



#[export]
lud <- function(x) {
	.Call(Rfast2_lud, x)
}