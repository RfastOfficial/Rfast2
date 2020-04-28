#[export]
is.upper.tri <- function(x,diag = FALSE) {
	.Call(Rfast2_is_upper_tri,x,diag)
}