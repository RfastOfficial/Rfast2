boot.ttest1 <- function(x, m, R = 999) {
    
  n <- length(x)
  xbar <- sum(x)/n
  s <- Rfast::Var(x, std = TRUE)
  stat <- (xbar - m)/s
  z <- x - xbar + m
  zb <- matrix(sample(z, n * R, replace = TRUE), ncol = R)
  xbar <- Rfast::colmeans(zb)
  s <- Rfast::colVars(zb, std = TRUE)
  bstat <- (xbar - m)/s
  pvalue <- ( sum( abs(bstat) >= abs(stat) ) + 1 ) / (R + 1)
  res <- c(sqrt(n) * stat, pvalue)
  names(res) <- c("stat", "bootstrap p-value")
  res
}
  
