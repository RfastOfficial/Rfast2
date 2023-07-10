#[export]
trim.mean <- function(x, a = 0.05, parallel = FALSE) {
  .Call(Rfast2_trimmean, x, a, parallel)
}

#[export]
rowTrimMean <- function(x, a = 0.05, parallel = FALSE, cores = 0) {
  .Call(Rfast2_rowTrimMean, x, a, parallel, cores)
}

#[export]
colTrimMean <- function(x, a = 0.05, parallel = FALSE, cores = 0) {
  UseMethod("colTrimMean")
}

#[export s3]
colTrimMean.matrix <- function(x, a = 0.05, parallel = FALSE, cores = 0) {
  .Call(Rfast2_colTrimMean, x, a, parallel, cores)
}

#[export s3]
colTrimMean.data.frame <- function(x, a = 0.05, parallel = FALSE, cores = 0) {
  .Call(Rfast2_colTrimMean, x, a, parallel, cores)
}
