#[export]
simplex.mle <- function(x, tol = 1e-07) {
  n <- length(x)
  xx <- x * (1 - x)
  simplexfun <- function(m, xx, x)  sum( (x - m)^2 /xx ) / ( m^2 * (1 - m)^2 )
  mod <- optimise(simplexfun, c(0, 1), xx = xx, x = x, tol = tol)
  s <- sqrt( mod$objective/n )
  param <- c( mod$minimum, s)
  names(param) <- c("mean", "sigma")
  list(param = param, loglik = -0.5 * n * log(2 * pi) - 1.5 * sum( log(xx) ) - n * log(s) - n/2 )
}


# simplex.mle <- function (x, tol = 1e-09) {
#   n <- length(x)
#   xx <- x * (1 - x)
#   a <- min(x)
#   b <- max(x)
#   ratio <- 2/(sqrt(5) + 1)
#   m1 <- b - ratio * (b - a)
#   m2 <- a + ratio * (b - a)
#   f1 <-  - sum( (x - m1)^2 /xx ) / ( m1^2 * (1 - m1)^2 )
#   f2 <-  - sum( (x - m2)^2 /xx ) / ( m2^2 * (1 - m2)^2 )
#
#   while ( abs(f2 - f1) > tol ) {
#     if (f2 < f1) {
#       b <- m2
#       m2 <- m1
#       f2 <- f1
#       m1 <- b - ratio * (b - a)
#       f1 <-  - sum( (x - m1)^2 /xx ) / ( m1^2 * (1 - m1)^2 )
#     } else {
#       a <- m1
#       m1 <- m2
#       f1 <- f2
#       m2 <- a + ratio * (b - a)
#       f2 <-  - sum( (x - m2)^2 /xx ) / ( m2^2 * (1 - m2)^2 )
#     }
#   }
#
#   m <- 0.5 * (m1 + m2)
#   s <- sqrt(  - f2/n )
#   param <- c( m, s)
#   names(param) <- c("mean", "sigma")
#   list(param = param, loglik = -0.5 * n * log(2 * pi) - 1.5 * sum( log(xx) ) - n * log(s) - n/2 )
# }
