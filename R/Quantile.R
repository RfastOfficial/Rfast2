#[export]
Quantile <- function(x, probs) {
	.Call(Rfast2_Quantile, x, probs)
}

#[export]
rowQuantile <- function(x, probs, parallel = FALSE) {
	.Call(Rfast2_rowQuantile, x, probs, parallel)
}

#[export s3]
colQuantile.data.frame<-function(x, probs, parallel = FALSE) {
	.Call(Rfast2_colQuantile, x, probs, parallel)
}
#[export]
colQuantile<-function(x, probs, parallel = FALSE) {
	UseMethod("colQuantile")
}
#[export s3]
colQuantile.matrix<-function(x, probs, parallel = FALSE) {
	.Call(Rfast2_colQuantile, x, probs, parallel)
}
