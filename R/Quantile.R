#[export]
Quantile <- function(x, probs) {
  .Call(Rfast2_Quantile, x, probs)
}

#[export]
rowQuantile <- function(x, probs, parallel = FALSE, cores = 0) {
  .Call(Rfast2_rowQuantile, x, probs, parallel, cores)
}

#[export s3]
colQuantile.data.frame <- function(x, probs, parallel = FALSE, cores = 0) {
  .Call(Rfast2_colQuantile, x, probs, parallel, cores)
}
#[export]
colQuantile <- function(x, probs, parallel = FALSE, cores = 0) {
  UseMethod("colQuantile")
}
#[export s3]
colQuantile.matrix <- function(x, probs, parallel = FALSE, cores = 0) {
  .Call(Rfast2_colQuantile, x, probs, parallel, cores)
}
