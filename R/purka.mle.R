#[export]
purka.mle <- function(x, tol = 1e-07) {
  if ( !is.matrix(x) )  x <- cbind( cos(x), sin(x) )
  p <- dim(x)[2]

  theta <- Rfast::mediandir(x)
  a <- x %*% theta
  a[ abs(a) > 1 ] <- 1
  A <- sum( acos(a) )
  n <- dim(x)[1]
  circle <- function(a, A, n)  n * log(a) - n * log(2) - n * log( 1 - exp( - a * pi ) ) - a * A
  sphere <- function(a, A, n)  n * log(a^2 + 1) - n * log(2 * pi) - n * log( 1 + exp( - a * pi ) ) - a * A
  hypersphere <- function(a, A, n) {
    n * lgamma(p/2) - 0.5 * n * p * log(pi) + n * ( log(besselI(a, p - 1, expon.scaled = TRUE)) + a ) - a * A
  }

  if (p == 2) {
    lika <- optimize(circle, c(0.001, 30000), maximum = TRUE, A = A, n = n, tol = tol) 
    a <- lika$maximum  ## estimated kappa
    f2 <-  -n / a^2 + n * pi^2 * exp( -a * pi)/ ( 1 - exp(-a * pi) )^2
  } else if (p == 3) {
    lika <- optimize(sphere, c(0.001, 30000), maximum = TRUE, A = A, n = n, tol = tol) 
    a <- lika$maximum  ## estimated kappa
    f2 <-  - (2 * a^2 * n - 2 * n) / (a^2 + 1)^2 - n * pi^2 * ( 1 + exp( a * pi) )^(-2) * exp( a * pi )    
  } else {
    lika <- optimize(hypersphere, c(0.001, 30000), maximum = TRUE, A = A, n = n, tol = tol)
    a <- lika$maximum  ## estimated kappa
    up1 <- ( besselI(a, p - 3) + 2 * besselI(p - 1, a) + besselI(p + 1, a) ) * besselI(p - 1, a)
    up2 <- ( besselI(a, p - 2) + 2 * besselI(p, a) )^2
    f2 <- 0.5 * n * ( up1 - up2 ) / besselI(p - 1, a)^2
  }
  ## f2 is the second derivative of the log-likelihood w.r.t alpha
    
  list( theta = theta, alpha = a, loglik = lika$objective, alpha.sd = 1 / sqrt( - f2 ) ) 
}
  
  
  

