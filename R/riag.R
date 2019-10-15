#[export]
riag <- function(n, mu) {
  p <- length(mu)
  x <- Rfast::matrnorm(n, p) +  rep(mu, rep(n, p))
  x / sqrt( Rfast::rowsums(x^2) )
}
