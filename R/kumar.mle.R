#[export]
kumar.mle <- function(x, tol = 1e-07, maxiters = 50) {

  n <- length(x)
  lx <- log(x)
  slx <- sum(lx) 
  
  ini <- Rfast::beta.mle(x)$param
  expa <- ini[1]   ;   expb <- ini[2]
  xa <- x^expa
  ya <- 1 - xa
  com <- xa * lx / ya
  scom <- sum(com)
  derab <-  - expb * expa * scom
  dera <- n + expa * slx + (1 - 1/ expb) * derab 
  dera2 <- expa * slx - (expb - 1) * expa^2 * sum( com * lx / ya )
  derb2 <- expb * sum( log(ya) )
  derb <- n + derb2
  aold <- c( log(expa), log(expb) ) 
  anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  i <- 2
  while ( sum( abs(anew - aold) ) > tol  &  i < maxiters ) {
    i <- i + 1
    aold <- anew
    expa <- exp( aold[1] )     ;      expb <- exp( aold[2] )
    xa <- x^expa
    ya <- 1 - xa
    com <- xa * lx / ya
    scom <- sum(com)
    derab <-  - expb * expa * scom
    dera <- n + expa * slx + (1 - 1/ expb) * derab  
    dera2 <- expa * slx - (expb - 1) * expa^2 * sum( com * lx / ya )
    derb2 <- expb * sum( log(ya) )
    derb <- n + derb2
    anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  }

  a <- exp( anew[1] )    ;     b <- exp( anew[2] )
  param <- c(a, b)
  loglik <-  n * log(a * b) + (a - 1) * slx + (b - 1) * derb2/expb
  names(param) <- c("shape", "scale")
  list(iters = i, param = param, loglik = loglik) 
}



#kumar <- function(x) {
#  n <- length(x)
#  slx <- sum( log(x) )
#  fa <- function(pa) {
#     a <- exp( pa[1] )    ;    b <- exp( pa[2] )
#     - n * log(a * b) - (a - 1) * slx - (b - 1) * sum( log(1 - x^a) )
#  }
#  ini <- log( Rfast::beta.mle(x)$param )
#  options(warn = -1)
#  mod <- nlm( fa, ini )
#  list( param = exp(mod$estimate), loglik = - mod$minimum )
#}





