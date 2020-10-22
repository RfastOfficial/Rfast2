mci <- function(fun, R = 10^6) {
  x <- Rfast::Rnorm(R)
  mean( fun(x) )
}

 