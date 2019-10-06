#[export]
censpois.mle <- function(x, tol = 1e-07) {
  ca <- min(x)
  z <- 0:ca
  fac <- factorial(z)
  n <- length(x)
  n2 <- sum(x == ca)
  n1 <- n - n2
  sx <- sum( x[x > ca] )
  a1 <- log( sx/n )
  down <- sum( exp(a1 * z)/fac )
  dera <-  - n1 * exp(a1) + sx + n2 * sum( z * exp( a1 * z )/fac ) / down
  dera2 <-  - n1 * exp(a1) + 
        n2 * ( sum( z^2 * exp( a1 * z/fac ) ) - sum( z * exp(a1 * z)/fac )^2 ) / down^2
  a2 <- a1 - dera/dera2
  i <- 2
  while ( abs(a2 - a1) > tol ) {
    i <- i + 1
    down <- sum( exp(a1 * z)/fac )
    a1 <- a2
    dera <-  - n * exp(a1) + sx + n2 * sum( z * exp( a1 * z )/fac ) / down
    dera2 <-  - n * exp(a1) + 
        n2 * ( sum( z^2 * exp( a1 * z/fac ) ) - sum( z * exp(a1 * z)/fac )^2 ) / down^2
    a2 <- a1 - dera/dera2  
  }

  loglik <-  - n * exp(a1) + sx * a1 - sum( lgamma(x + 1) ) + n2 * log( down )
  list(iters = i, loglik = loglik, lambda = exp(a2) )
}
