#[export]
moment <- function(x, order = 1, central = FALSE, absolute = FALSE)
  n <- length(x)
  if ( central )   x <- x - sum(x)/n
  if ( absolute )  x <- abs(x)
  sum( x^order ) / n
}


#[export]
moments <- function(x, order = 1, central = FALSE, absolute = FALSE)
  n <- dim(x)[1]
  if ( central )   x <- Rfast::eachrow(x, Rfast::colmeans(x), oper = "-")
  if ( absolute )  x <- abs(x)
  Rfast::colmeans( x^order ) 
}