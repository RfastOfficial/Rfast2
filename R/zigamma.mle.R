#[export]
zigamma.mle <- function(x, tol = 1e-07) {
  n <- length(x)
  x1 <- x[x > 0]
  n1 <- length(x1)
  n0 <- n - n1
  prob <- n1/n
  lik0 <- n0 * log(1 - prob) + n1 * log(prob)
  mod <- Rfast::gammamle(x1, tol = tol)
  param <- c(prob, mod$param)
  names(param) <- c("prop1", "shape", "scale")
  list(iters = mod$iters, loglik = sum(lik0, mod$loglik, na.rm = TRUE), param = param)
}
