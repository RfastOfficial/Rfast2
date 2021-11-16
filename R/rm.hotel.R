rm.hotel <- function(x, a = 0.05) {
  ## x is the data set
  ## a is the level of significance set by default to 0.05
  m <- Rfast::colmeans(x)
  s <- cov(x)  ## sample mean vector and covariance matrix
  n <- dim(x)[1]  ## sample size
  p <- dim(x)[2]  ##  dimensionality of the data
  C <-  - diag(p)
  C[, 1] <- 1
  A <- C %*% m
  B <- solve( (C %*% s) %*% C, A)
  T2 <- n * sum(A * B)
  test <- (n - p + 1) / (n - 1) / (p - 1) * T2  ## test statistic
  pvalue <- pf(test, p - 1, n - p + 1, lower.tail = FALSE)  ## p-value of the test statistic
  crit <- qf(1 - a, p - 1, n - p + 1)  ## critical value of the F disitribution
  result <- c(test, pvalue, crit, p - 1, n - p + 1)
  names(result) <- c("test", "p-value", "numer df", "denom df", "critical")
  list(m = m, result = result)
}
