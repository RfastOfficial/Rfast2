#[export]
is.lower.tri <- function(x,diag = FALSE) {
	.Call(Rfast2_is_lower_tri,x,diag)
}