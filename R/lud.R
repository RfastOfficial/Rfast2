#[export]
lud <- function(x) {
	.Call(Rfast2_lud,x)
}