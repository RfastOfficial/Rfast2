#[export]
zil.mle <- function(x) {

  n <- length(x)
  x1 <- x[x > 0]
  n1 <- length(x1)
  n0 <- n - n1
  prob <- n1/n
  lik0 <- n0 * log(1 - prob) + n1 * log(prob)

  lx1 <- log(x1)
  lx2 <- log(1 - x1)
  y <- lx1 - lx2
  sy <- sum(y)
  m <- sy/n1
  s <- (sum(y^2) - n1 * m^2)/n1
  loglik <- sum(dnorm(y, m, sqrt(s), log = TRUE)) - sy
  param <- c(prob, m, n1 * s/(n1 - 1))
  names(param) <- c("prop1", "mean", "unbiased variance")
  list(loglik = sum(lik0, loglik, na.rm = TRUE), param = param)
}
