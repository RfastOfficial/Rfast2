#[export]
walter.ci <- function(x1, x2, n1, n2, a = 0.05) {
  x <- x1 + 0.5
  y <- x2 + 0.5
  n <- n1 + 0.5
  m <- n2 + 0.5
  se <- sqrt( 1/x + 1/y - 1/m - 1/n )
  lr <- log(x/n) - log(y/m)
  ci <- c(lr - qnorm(1 - a/2) * se, lr + qnorm(1 - a/2) * se)
  list( rat = exp(lr), ci = exp(ci) )
} 




