#[export]
collogitnorm.mle <- function(x) {
  n <- dim(x)[1]
  lx1 <- Rfast::Log(x)
  lx2 <- Rfast::Log(1 - x)
  y <- lx1 - lx2
  sy <- Rfast::colsums(y)
  m <- sy/n
  s <- ( Rfast::colsums(y^2) - n * m^2 ) / n
  loglik <- Rfast::rowsums( dnorm(t(y), m, sqrt(s), log = TRUE) ) - sy
  res <- cbind(m, n * s/(n - 1), loglik)
  colnames(res) <- c("mean", "unbiased variance", "loglik")
  res
}
