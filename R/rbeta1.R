#[export]
rbeta1 <- function(n, a) {
  if (a != 1 ) {
    x <- exp( -rexp(n, a) )
  } else x <- Runif(n)
  x   
}