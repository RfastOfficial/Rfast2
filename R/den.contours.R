#[export]
den.contours <- function(x, type = "normal", v = 5) {
  n1 <- 100 
  n2 <- 100
  x1 <- seq( min(x[, 1]) - 0.5, max(x[, 1]) + 0.5, length = n1 )
  x2 <- seq( min(x[, 2]) - 0.5, max(x[, 2]) + 0.5, length = n2 )
  ## if for example the y axis is longer than the x axis, then you might 
  ## want to change n2

  if ( type == "normal" ) {
    m <- Rfast::colmeans(x) ## mean vector  
    s <- cov(x) ## covariance matrix
    r <- s[2] / sqrt(s[1] * s[4])
    con <-  - log(2 * pi) - log( det(s) )  ## constant part
    z1 <- ( x1 - m[1] ) / sqrt( s[1] )
    z2 <- ( x2 - m[2] ) / sqrt( s[4] )
    mat1 <- outer(z1^2, rep(1, n1), "*")
    mat2 <- outer(rep(1, n2), z2^2,  "*")
    mat3 <- tcrossprod(z1, z2)
    mat <- con - 0.5 /(1 - r^2) * (mat1 + mat2 - 2 * r * mat3)
   	
  } else  if ( type == "t" ) {
    ## we will use the previous function 'multivt' to 
    ## estimate the parameters of the bivariate t first
    f <- Rfast::mvt.mle(x) 
    m <- f$location
    s <- f$scatter
    r <- s[2] / sqrt(s[1] * s[4])
    con <- lgamma( (v + 2) / 2 ) - lgamma(v / 2) - 0.5 * log( det(pi * v * s) )
    z1 <- ( x1 - m[1] ) / sqrt(s[1])
    z2 <- ( x2 - m[2] ) / sqrt(s[4])
    mat1 <- outer(z1^2, rep(1, n1), "*")
    mat2 <- outer(rep(1, n2), z2^2,  "*")
    mat3 <- tcrossprod(z1, z2)
    mat <- con - 0.5 * (v + 2) * log1p( 1 /(1 - r^2) * 
    (mat1 + mat2 - 2 * r * mat3) / v )
   
  } else if ( type == "lnorm" ) {
    y <- Rfast::Log(x)
    m <- Rfast::colmeans(y)  ## fitted mean vector of the logged data
    s <- cov(y) ## estimated covariance matrix of the logged data 
    r <- s[2] / sqrt(s[1] * s[4])
    con <-  -log(2 * pi) - 0.5 * log(det(s))  ## contant part
    x1 <- seq( max(min(x[, 1]) - 0.5, 0.01), max(x[, 1]) + 0.5, length = n1 )
    x2 <- seq( max(min(x[, 2]) - 0.5, 0.01), max(x[, 2]) + 0.5, length = n1 )
    z1 <- ( log(x1) - m[1] ) / sqrt(s[1])
    z2 <- ( log(x2) - m[2] ) / sqrt(s[4])
    xat1 <- matrix( rep(x1, n2), ncol = n2 )
    xat2 <- matrix( rep(x2, n1), ncol = n2, byrow = TRUE )
    mat3 <- tcrossprod(z1, z2)
    mat <- con - xat1 - xat2 - 0.5 /(1 - r^2) * (mat1 + mat2 - 2 * r * mat3)
    m <- exp( m + 0.5 * diag(s) )
  }    

  mat <- exp(mat)
  ind <- ( mat < Inf )  
  ind[ ind == FALSE ] <- NA
  mat <- mat * ind   
  ## we did this to avoid any issues with high numbers
  contour(x1, x2, mat, nlevels = 10, col = 2, xlab = colnames(x)[1], 
  ylab = colnames(x)[2], cex.lab = 1.2, cex.axis = 1.2)
  points(x[, 1], x[, 2], pch = 20)
  points(m[1], m[2], pch = 10, cex = 1.5)
}