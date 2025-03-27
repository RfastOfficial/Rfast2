
#[export]
kernel <- function(x, h = "silverman", parallel = FALSE, cores = 0) {
  res <- .Call(Rfast2_kernel,t(x),h, parallel, cores)
  if(is.matrix(res)) res <- t(res)
  res
}


.kernel <- function(x, h = "silverman") {
  
  if ( is.vector(x) ) {
    n <- length(x)
    lenh <- length(h)

    if ( lenh == 1 ) {
      if ( h == "silverman" ) {
        s <- Rfast::Var(x, std = TRUE)
        iqr <- diff( Rfast::Quantile(x, probs = 0.25, 0.75) )
        h <- 0.9 * min(s, iqr / 1.34 ) * n^( -0.2 )
      } else if ( h == "scott" ) {
        s <- Rfast::Var(x, std = TRUE)
        h <- 1.06 * s * n^( -0.2 )
      }
      h2 <- 2 * h^2
      d <- Rfast::Dist(x, square = TRUE) / h2
      f <- ( Rfast::colsums( exp(-d) ) - 1) / ( (n - 1) * h * sqrt(2 * pi) )

    } else {  ## h is a vector
      f <- matrix( nrow = n, ncol = lenh )         
      h2 <- 2 * h^2
      d <- Rfast::Dist(x, square = TRUE) 
      for ( j in 1:lenh ) {
        f[, j] <- ( Rfast::colsums( exp(-d / h2[j]) ) - 1) / ( (n - 1) * h[j] * sqrt(2 * pi) )
      } 
    }

  } else { ## x is a matrix
    n <- dim(x)[1]  ;  p <- dim(x)[2]

    if ( h == "silverman" ) {
      s <- Rfast::colVars(x, std = TRUE)
      iqr <- Rfast::colQuantile(x, probs = 0.25, 0.75)
      iqr <- iqr[2, ] - iqr[1, ]
      h <- 0.9 * pmin(s, iqr / 1.34 ) * n^( -0.2 )
    } else if ( h == "scott" ) {
      s <- Rfast::colVars(x, std = TRUE)
      h <- 1.06 * s * n^( -0.2 )
    }
      
    x <- t( t(x) / ( sqrt(2) * h ) )
    d <- Rfast::Dist(x, square = TRUE)
    f <- ( Rfast::colsums( exp(-d) ) - 1) / ( (n - 1) * prod(h) * (2 * pi)^(0.5 * p) )
  }
 
  f
} 



