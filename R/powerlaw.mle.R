#[export]
powerlaw.mle <- function(x) {
  n <- length(x)
  x1 <- min(x)
  com <- sum( log(x) ) - n * log(x1)
  a <- 1 + n / com
  loglik <- n * log( (a - 1) / x1 ) - a * com 
  list(alpha = a, loglik = loglik)
}