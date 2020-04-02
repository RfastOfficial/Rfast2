#[export]
js.propreg <- function(y, x, tol = 1e-07, maxiters = 100) {
  
  x <- model.matrix(~., data.frame(x) )
  n <- dim(x)[1]
  p <- mean(y)
  be <- c( log(p / (1 - p ) ), numeric(dim(x)[2] - 1) )
  d1 <- sum( (y + p) * log(y + p) ) + n * p * log(2 * p) + n * (1 - p) * log(1 - p) - sum( (2 - y - p) * log(2 - y - p) )
 
  .e3 <-  1/p - 1   ## as.vector( exp(-x %*% be) )
  .e4 <- 1 + .e3
  .e5 <- .e4^2
  .e6 <- 1/.e4
  .e7 <- .e6 + y
  com <- ( - log(1 - .e6) + log(.e7) + 2.693147 - log1p(.e3) + log(2 - .e7) ) * .e3/.e5 
  der <- Rfast::eachcol.apply(x, com) 

  ## der <-  -Rfast::eachcol.apply( x, (1 + log(1 - .e6)) * com ) + Rfast::eachcol.apply( x, (1 + log(.e7)) * com ) + 
  ## Rfast::eachcol.apply( x, (1.693147 - log1p(.e3)) * com) + Rfast::eachcol.apply( x, ( 1 + log(2 - .e7) ) * com )
  
  .e5 <- 1/.e4
  .e6 <- .e5 + y
  .e7 <- .e4^2
  .e9 <- 1 - .e5
  .e10 <- 2 - .e6
  .e11 <- 1 + log(.e6)
  .e12 <- log(.e9)
  .e13 <- log(.e10)
  .e14 <- log1p(.e3)
  com34 <- .e3/.e4
  com37 <- .e3/.e7
  com <- .e11 + .e13 -.e12 - .e14
  a2 <- ( ( 4.386294 + 2 * com + 1/(.e9 * .e4) + 1/(.e4 * .e6) - 1/(.e4 * .e10) ) * com34 - 1.693147 - com ) * com37
  der2 <- crossprod(x * a2, x)

  ## der2 <- crossprod( x * ( (1 + 2 * (1.693147 - .e14)) * com34 + .e14 - 1.693147 ) * com37, x) + 
  ## crossprod( x * ( (1/(.e9 * .e4) - 2 * (1 + .e12)) * com34 + 1 + .e12) * com37, x ) + 
  ## crossprod( x * ( (1/(.e4 * .e6) + 2 * .e11) * com34 - .e11) * com37, x ) - 
  ## crossprod( x * ( (1/(.e4 * .e10) - 2 * (1 + .e13)) * com34 + 1 + .e13) * com37, x )

  be <- be - solve(der2, der)
  p <- 1 / ( 1 + exp(-x %*% be) )
  d2 <- sum( (y + p) * log(y + p) ) + sum( p * log(2 * p) ) + sum( (1 - p) * log(1 - p) ) - 
  sum( (2 - y - p) * log(2 - y - p) )

  i <- 2
  while ( d1 - d2 > tol  &  i < maxiters) {
    i <- i + 1
    d1 <- d2
    .e3 <- as.vector( exp(-x %*% be) )
    .e4 <- 1 + .e3
    .e5 <- .e4^2
    .e6 <- 1/.e4
    .e7 <- .e6 + y
    com <- ( - log(1 - .e6) + log(.e7) + 2.693147 - log1p(.e3) + log(2 - .e7) ) * .e3/.e5 
    der <- Rfast::eachcol.apply(x, com) 

    .e5 <- 1/.e4
    .e6 <- .e5 + y
    .e7 <- .e4^2
    .e9 <- 1 - .e5
    .e10 <- 2 - .e6
    .e11 <- 1 + log(.e6)
    .e12 <- log(.e9)
    .e13 <- log(.e10)
    .e14 <- log1p(.e3)
    com34 <- .e3/.e4
    com37 <- .e3/.e7
    com <- .e11 + .e13 -.e12 - .e14
    a2 <- ( ( 4.386294 + 2 * com + 1/(.e9 * .e4) + 1/(.e4 * .e6) - 1/(.e4 * .e10) ) * com34 - 1.693147 - com ) * com37
    der2 <- crossprod(x * a2, x)
   
    be <- be - solve(der2, der)
    p <- 1 / ( 1 + exp(-x %*% be) )
    d2 <- sum( (y + p) * log(y + p) ) + sum( p * log(2 * p) ) + sum( (1 - p) * log(1 - p) ) - 
    sum( (2 - y - p) * log(2 - y - p) )

  }
  list(be = be, der2 = der2, js = d2, iters = i) 
}




#[export]
helinger.propreg(y, x, tol = 1e-07, maxiters = 100) {

  x <- model.matrix( y~., data.frame(x) )
  sqy <- sqrt(y)
  dm <- dim(x)
  n <- dm[1]    ;    d <- dm[2]
  p <- sum(y)/n
  be <- c( log( p / (1 - p) ), numeric(d - 1) )
  d1 <- sum( ( sqy - sqrt(p) )^2 )
  
  .e3 <- as.vector( exp(-x %*% be) )
  .e4 <- 1 + .e3
  .e6 <- sqrt(1/.e4)
  der <-  - Rfast::eachcol.apply(x, .e3 * (sqy - .e6) / (.e4^2 * .e6) )

  .e8 <- .e4^2 * .e6
   der2 <- crossprod(x * ( ( 0.5 * (.e3/.e8) + sqy - .e6 )/.e8 + ( 0.5/(.e4 * .e6) - 2 * .e6 ) * 
           .e4 * .e3 * ( sqy - .e6)/.e8^2 ) * .e3, x)
   be <- be - solve(der2, der)
   .e3 <- as.vector( exp(-x %*% be) )
   .e4 <- 1 + .e3
   .e6 <- sqrt(1/.e4)
   d2 <- sum( ( sqy - .e6 )^2 )
   i <- 2

  while ( d1 - d2 > tol  &  i < maxiters ) {
    i <- i + 1
    d1 <- d2
    der <-  - Rfast::eachcol.apply(x, .e3 * (sqy - .e6) / (.e4^2 * .e6) )
    .e8 <- .e4^2 * .e6
    der2 <- crossprod(x * ( ( 0.5 * (.e3/.e8) + sqy - .e6 )/.e8 + ( 0.5/(.e4 * .e6) - 2 * .e6 ) * 
           .e4 * .e3 * ( sqy - .e6)/.e8^2 ) * .e3, x)
    be <- be - solve(der2, der)
    .e3 <- as.vector( exp(-x %*% be) )
    .e4 <- 1 + .e3
    .e6 <- sqrt(1/.e4)
    d2 <- sum( ( sqy - .e6 )^2 )
  } 

  list(be = be, der2 = der2, H = d2, iters = i) 
}

  


     
