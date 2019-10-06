
#[export]
Intersect <- function(x,y) {
	.Call(Rfast2_inter,x,y)
}