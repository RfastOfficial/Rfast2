#[export]
halfcauchy.mle <- function(x, tol = 1e-07, maxiters = 50) {
  n <- length(x)
  ea <- 0.5 * ( Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4) )
  a <- log(ea)
  x2 <- x^2
  com <- 1 / (x2 + ea^2 )
  lik1 <- n * a + sum( log(com) )
  der <- n - 2 * ea^2 * sum(com)  
  der2 <-  - 4 * ea^2 * sum(x2 * com^2)
  a <- a - der/der2
  ea <- exp(a)
  com <- 1 / (x2 + ea^2 )
  lik2 <- n * a + sum( log(com) )
  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters ) {
    i <- i + 1
    lik1 <- lik2
    der <- n - 2 * ea^2 * sum(com)  
    der2 <-  - 4 * ea^2 * sum(x2 * com^2)
    a <- a - der/der2
    ea <- exp(a)
    com <- 1 / (x2 + ea^2 )
    lik2 <- n * a + sum( log(com) )
  }  
  list(iters = i, loglik = lik2 + n * log(2/pi), scale = ea)
}

  

