#[export]
Merge <- function(x,y) {
	.Call(Rfast2_merge,x,y)
}