#[export]
pinar1 <- function(x, unbiased = FALSE) {
  n <- length(x)
  n1 <- n - 1
  sumx1 <- sum(x[-1])
  sumxn <- sum(x[-n]) 
  num <- sum( x[-1] * x[-n] ) - sumx1 / n1 * sumxn
  denom <- sum(x[-n]^2) - sumxn / n1 * sumxn
  alpha <- num / denom
  if ( unbiased )  alpha <-  ( n * alpha + 1 ) / (n - 3) 
  lambda <- ( sumx1 - alpha * sumxn ) / n1
  res <- c(lambda, alpha)
  names(res) <- c("lambda", "alpha")
  res
}


#[export]
colpinar1 <- function(x, unbiased = FALSE) {
  n <- dim(x)[1]
  n1 <- n - 1
  sumx1 <- Rfast::colsums(x[-1, ])
  sumxn <- Rfast::colsums(x[-n, ]) 
  num <- Rfast::colsums( x[-1, ] * x[ -n, ] ) - sumx1 / n1 * sumxn
  denom <- Rfast::colsums(x[-n, ]^2 ) - sumxn / n1 * sumxn
  alpha <- num / denom
  if ( unbiased )  alpha <-  ( n * alpha + 1 ) / (n - 3) 
  lambda <- ( sumx1 - alpha * sumxn ) / n1
  res <- cbind(lambda, alpha)
  colnames(res) <- c("lambda", "alpha")
  res
}
