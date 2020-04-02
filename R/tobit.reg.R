#[export]
tobit.reg <- function(y, x, ylow = 0, full = FALSE, tol = 1e-07, maxiters = 100) {

  x <- model.matrix(y ~., data.frame(x) )
  ind <- which( y == ylow)
  y1 <- y[ -ind ]
  y0 <- y[ ind ]
  x1 <- x[ -ind, ]   
  x0 <- x[ ind, ]
  n <- length(y)
  n1 <- length(y1)
   
  be <- solve( crossprod(x), crossprod(x, y) )
  e1 <- y1 - as.vector( x1 %*% be )
  e0 <- y0 - as.vector( x0 %*% be )
  a <- 0.5 * log( ( sum(e1^2) + sum(e0^2) ) / n )  

  s <- exp(a)
  lik1 <-  - n1 * log(s) + sum( log( dnorm( e1 / s) ) ) + sum( log( pnorm( e0 / s) ) )   

  .e1 <- exp(a)
  .e3 <- as.vector( y1 - x1 %*% be )
  .e4 <- as.vector(y0 - x0 %*% be)/.e1
  .e5 <- .e3/.e1
  derb <-  - Rfast::eachcol.apply( x0, dnorm(.e4, 0, 1) / ( .e1 * pnorm(.e4) ) ) + 
             Rfast::eachcol.apply( x1, dnorm(.e5, 0, 1) * .e3 / ( dnorm(.e5) * .e1^2 ) )
  .e4 <- .e4 * .e1 
  .e6 <- .e4 / .e1
  .e7 <- .e3 / .e1
  dera <-  - sum( dnorm(.e6, 0, 1) * .e4/(.e1 * pnorm(.e6)) ) + sum( dnorm(.e7, 0, 1) * .e3^2/(dnorm(.e7) * .e1^2) ) - n1
  
  .e5 <- .e4 
  .e6 <- .e5/.e1
  .e8 <- .e1^2
  .e9 <- .e3^2
  .e10 <- dnorm(.e6, 0, 1)
  .e12 <- dnorm(.e7) * .e8
  .e13 <- dnorm(.e7, 0, 1)
  .e14 <- pnorm(.e6)
  derb2 <- -crossprod( x0 * (.e5/(.e8 * .e14) + .e10 * .e1/(.e1 * .e14)^2) * .e10/.e1, x0) + 
            crossprod(x1 * ( (.e9/.e8 - 1)/.e12 - .e13 * .e9/.e12^2 ) * .e13, x1)

  .e4 <- y0 - x0 %*% be
  .e5 <- y1 - x1 %*% be
  .e6 <- .e4/.e1
  .e7 <- .e5/.e1
  .e8 <- .e5^2
  .e9 <- dnorm(.e7)
  .e10 <- pnorm(.e6)
  .e11 <- dnorm(.e6, 0, 1)
  .e12 <- dnorm(.e7, 0, 1)
  .e13 <- .e1 * .e10
  dera2 <-  - sum( ( ( (2 * (.e9 * .e1) + .e12 * .e8/.e1) * .e1/(.e9 * .e1^2)^2 - .e8/(.e9 * .e1^4) ) * .e12 * .e8 ) ) + 
              sum( ( (.e13 - .e11 * .e4)/.e13^2 - .e4^2/(.e1^3 * .e10) ) * .e11 * .e4 )

  .e3 <- as.vector( y1 - x1 %*% be )
  .e5 <- as.vector( y0 - x0 %*% be )
  .e6 <- .e5/.e1
  .e7 <- .e3/.e1
  .e8 <- dnorm(.e7)
  .e9 <- pnorm(.e6)
  .e10 <- .e3^2
  .e11 <- dnorm(.e6, 0, 1)
  .e12 <- dnorm(.e7, 0, 1)
  .e13 <- .e1 * .e9
  derab <-  - colsums( x1 * ( (2 * (.e8 * .e1) + .e12 * .e10/.e1 ) * .e1 / ( .e8 * .e1^2)^2 - .e10 / (.e8 * .e1^4) ) * .e12 * .e3 ) + 
  colsums( x0 * ( (.e13 - .e11 * .e5)/.e13^2 - .e5^2/(.e1^3 * .e9) ) * .e11 )
  
  der <- c(dera, derb)
  der2 <- cbind(derab, derb2)
  der2 <- rbind( c(dera2, derab), der2 )
  
  par <- c(a, be) - solve(der2, der)
  a <- par[1]  ;  be <- par[-1]
  e1 <- y1 - as.vector( x1 %*% be )
  e0 <- y0 - as.vector( x0 %*% be )
  s <- exp(a)
  lik2 <-  - n1 * log(s) + sum( log( dnorm( e1 / s) ) ) + sum( log( pnorm( e0 / s) ) ) 

  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters ) {  
    i <- i + 1
    lik1 <- lik2
    .e1 <- exp(a)
    .e3 <- as.vector( y1 - x1 %*% be )
    .e4 <- as.vector(y0 - x0 %*% be)/.e1
    .e5 <- .e3/.e1
    derb <-  - Rfast::eachcol.apply( x0, dnorm(.e4, 0, 1) / ( .e1 * pnorm(.e4) ) ) + 
               Rfast::eachcol.apply( x1, dnorm(.e5, 0, 1) * .e3 / ( dnorm(.e5) * .e1^2 ) )
    .e4 <- .e4 * .e1 
    .e6 <- .e4 / .e1
    .e7 <- .e3 / .e1
    dera <-  - sum( dnorm(.e6, 0, 1) * .e4/(.e1 * pnorm(.e6)) ) + sum( dnorm(.e7, 0, 1) * .e3^2/(dnorm(.e7) * .e1^2) ) - n1
  
    .e5 <- .e4 
    .e6 <- .e5/.e1
    .e8 <- .e1^2
    .e9 <- .e3^2
    .e10 <- dnorm(.e6, 0, 1)
    .e12 <- dnorm(.e7) * .e8
    .e13 <- dnorm(.e7, 0, 1)
    .e14 <- pnorm(.e6)
    derb2 <-  -crossprod( x0 * (.e5/(.e8 * .e14) + .e10 * .e1/(.e1 * .e14)^2) * .e10/.e1, x0) + 
               crossprod(x1 * ( (.e9/.e8 - 1)/.e12 - .e13 * .e9/.e12^2 ) * .e13, x1)

    .e4 <- y0 - x0 %*% be
    .e5 <- y1 - x1 %*% be
    .e6 <- .e4/.e1
    .e7 <- .e5/.e1
    .e8 <- .e5^2
    .e9 <- dnorm(.e7)
    .e10 <- pnorm(.e6)
    .e11 <- dnorm(.e6, 0, 1)
    .e12 <- dnorm(.e7, 0, 1)
    .e13 <- .e1 * .e10
    dera2 <-  - sum( ( ( (2 * (.e9 * .e1) + .e12 * .e8/.e1) * .e1/(.e9 * .e1^2)^2 - .e8/(.e9 * .e1^4) ) * .e12 * .e8 ) ) + 
                sum( ( (.e13 - .e11 * .e4)/.e13^2 - .e4^2/(.e1^3 * .e10) ) * .e11 * .e4 )

    .e3 <- as.vector( y1 - x1 %*% be )
    .e5 <- as.vector( y0 - x0 %*% be )
    .e6 <- .e5/.e1
    .e7 <- .e3/.e1
    .e8 <- dnorm(.e7)
    .e9 <- pnorm(.e6)
    .e10 <- .e3^2
    .e11 <- dnorm(.e6, 0, 1)
    .e12 <- dnorm(.e7, 0, 1)
    .e13 <- .e1 * .e9
    derab <-  - colsums( x1 * ( (2 * (.e8 * .e1) + .e12 * .e10/.e1 ) * .e1 / ( .e8 * .e1^2)^2 - .e10 / (.e8 * .e1^4) ) * .e12 * .e3 ) + 
    colsums( x0 * ( (.e13 - .e11 * .e5)/.e13^2 - .e5^2/(.e1^3 * .e9) ) * .e11 )
  
    der <- c(dera, derb)
    der2 <- cbind(derab, derb2)
    der2 <- rbind( c(dera2, derab), der2 )
  
    par <- c(a, be) - solve(der2, der)
    a <- par[1]  ;  be <- par[-1]
    e1 <- y1 - as.vector( x1 %*% be )
    e0 <- y0 - as.vector( x0 %*% be )
    s <- exp(a)
    lik2 <-  - n1 * log(s) + sum( log( dnorm( e1 / s) ) ) + sum( log( pnorm( e0 / s) ) ) 
  }
  names(lik2) <- NULL
  res <- list(be = be, s = s, loglik = lik2, iters = i)
  if ( full ) {
     se <- solve( der2 )
     se <- sqrt( diag(se)[-1] )
     wald <- be/se
     pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
     info <- cbind(be, se, wald, pval)
     colnames(info) <- c("estimate", "std error", "Wald stat", "p-value")
     rownames(info) <- colnames(x)
     res <- list(info = info, loglik = lik2)
  }
  res
}
   

  



