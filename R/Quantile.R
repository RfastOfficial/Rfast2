#[export]
Quantile <- function(x, probs) {
	.Call(Rfast2_Quantile, x, probs)
}




#[export]
rowQuantile <- function(x, probs, parallel = FALSE) {
	.Call(Rfast2_row_Quantile, x, probs, parallel)
}




#[export]
colQuantile<-function(x, probs, parallel = FALSE) {
	.Call(Rfast2_col_Quantile, x, probs, parallel)
}
