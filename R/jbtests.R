#[export]
jbtest <- function(x) {
  n <- length(x)
  S <- Rfast::skew(x)
  K <- Rfast::kurt(x)
  jb <- n / 6 * ( S^2 + 0.25 * (K - 3)^2 )
  pvalue <- pchisq(jb, 2, lower.tail = FALSE)
  res <- c(jb, pvalue)
  names(res) <- c("stat", "p-value")
  res
}


#[export]
jbtests <- function(x) {
  n <- dim(x)[1]
  S <- Rfast::colskewness(x)
  K <- Rfast::colkurtosis(x)
  jb <- n / 6 * ( S^2 + 0.25 * (K - 3)^2 )
  pvalue <- pchisq(jb, 2, lower.tail = FALSE)
  res <- cbind(jb, pvalue)
  colnames(res) <- c("stat", "p-value")
  res
}