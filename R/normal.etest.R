normal.etest <- function(x, R = 999) {
  n <- length(x)
  s <- Rfast::Var(x, std = TRUE)
  y <- ( x - mean(x) ) / s
  y <- sort(y)
  K <- seq(1 - n, n - 1, 2)
  stat <- 4 * sum( y * pnorm(y) + dnorm(y) ) - 2 * mean(K * y)  
  
  z <- Rfast::matrnorm(n, R)
  z <- Rfast::standardise(z)
  z <- Rfast::colSort(z)
  bootstat <- 4 * Rfast::colsums( z * pnorm(z) ) + 
              4 * Rfast::colsums( exp(-0.5 * z^2) ) / sqrt(2 * pi) -
              2 * Rfast::eachcol.apply(z, K) / n
  pvalue <- ( sum(bootstat > stat) + 1 ) / (R + 1)
  res <- c(stat - 2 * n / sqrt(pi), pvalue)
  names(res) <- c("statistic", "p-value")
  res
} 

  
 