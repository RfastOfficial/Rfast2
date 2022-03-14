
#[export]
Runif <- function(n, min = 0, max = 1){
	.Call(Rfast2_Runif, n, min, max)
}
#[export]
Sample <- function(x, size, replace = FALSE){
	.Call(Rfast2_Sample, x, size, replace)
}
#[export]
Sample.int <- function(n, size = n, replace = FALSE){
	.Call(Rfast2_Sample_int, n, size, replace)
}