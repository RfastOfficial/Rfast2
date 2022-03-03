#[export]
trim.mean <- function(x, a = 0.05) {
	.Call(Rfast2_trimmean, x, a)
}

#[export]
rowTrimMean <- function(x, a = 0.05, parallel = FALSE) {
	.Call(Rfast2_rowTrimMean, x, a, parallel)
}

#[export]
colTrimMean <- function(x, a = 0.05, parallel = FALSE) {
	UseMethod("colTrimMean")
}

#[export s3]
colTrimMean.matrix <- function(x, a = 0.05, parallel = FALSE) {
	.Call(Rfast2_colTrimMean, x, a, parallel)
}

#[export s3]
colTrimMean.data.frame <- function(x, a = 0.05, parallel = FALSE) {
	.Call(Rfast2_colTrimMean, x, a, parallel)
}
