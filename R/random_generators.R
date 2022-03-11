
#[export]
Runif <- function(n, min = 0, max = 1){
	.Call(Rfast2_Runif, n, min, max)
}
#[export]
Sample <- function(n, size = n, replace = FALSE){
	.Call(Rfast2_Sample, n, size, replace)
}
#[export]
Sample.int <- function(x, size, replace = FALSE){
	.Call(Rfast2_Sample_int, n, min, max)
}