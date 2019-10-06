#[export]
collognorm.mle <- function(x) {
  n <- dim(x)[1]
  x <- Rfast::Log(x)
  sx <- Rfast::colsums(x)
  m <- sx/n
  s <- Rfast::colsums(x^2)/n - m^2
  loglik <-  -0.5 * n * (log(2 * pi * s) + 1) - sx
  res <- cbind(m, s, loglik)
  colnames(res) <- c("mean", "variance", "loglik")
  res
}
