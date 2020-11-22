#[export]
censpois.mle <- function(x, tol = 1e-07) {
  ca <- min(x)
  z <- 0:ca
  fac <- factorial(z)
  n <- length(x)
  n2 <- sum(x == ca)
  n1 <- n - n2
  sx <- sum( x[x > ca] )
  a1 <- log( sx/n )
  down <- sum( exp(a1 * z)/fac )
  dera <-  - n1 * exp(a1) + sx + n2 * sum( z * exp( a1 * z )/fac ) / down
  dera2 <-  - n1 * exp(a1) + 
        n2 * ( sum( z^2 * exp( a1 * z/fac ) ) - sum( z * exp(a1 * z)/fac )^2 ) / down^2
  a2 <- a1 - dera/dera2
  i <- 2
  while ( abs(a2 - a1) > tol ) {
    i <- i + 1
    down <- sum( exp(a1 * z)/fac )
    a1 <- a2
    dera <-  - n * exp(a1) + sx + n2 * sum( z * exp( a1 * z )/fac ) / down
    dera2 <-  - n * exp(a1) + 
        n2 * ( sum( z^2 * exp( a1 * z/fac ) ) - sum( z * exp(a1 * z)/fac )^2 ) / down^2
    a2 <- a1 - dera/dera2  
  }

  loglik <-  - n * exp(a1) + sx * a1 - sum( lgamma(x + 1) ) + n2 * log( down )
  list(iters = i, loglik = loglik, lambda = exp(a2) )
}


#[export]
censweibull.mle <- function(x, di, tol = 1e-07) {
  y <- log(x)
  y1 <- y[di == 1]
  y2 <- y[di == 0]
  n <- length(x)
  n1 <- length(y1)
  n2 <- n - n1
  mod <- Rfast::weibull.mle(x[di==1])
  if ( n2 > 0 ) {
    m <- log(mod$param[2])
    es <- 1/mod$param[1]
    s <- log( es )
    z1 <- ( y1 - m ) / es
    z2 <- ( y2 - m ) /es
    com <- sum(z1)
    ez1 <- exp(z1) 
    ez2 <- exp(z2)
    com1 <- sum(ez1)
    com2 <- sum(ez2)
  
    lik1 <- com - com1 - n1 * s - com2

    derm2 <-  - com1 / es^2 - com2 / es^2  
    derm <-  - n1/ es + - derm2 * es 
    ders <-  - com + sum(ez1 * z1) - n1 + sum(ez2 * z2)
    ders2 <-  com - sum(ez1 * z1^2) - sum(ez1 * z1) - sum(ez2 * z2^2) - sum(ez2 * z2)
    derms <- n1/es - sum(ez1 * z1) / es - com1/es - sum(ez2 * z2) / es - com2/es 
    anew <- c(m, s) - c(ders2 * derm - derms * ders, -derms * 
              derm + derm2 * ders)/(derm2 * ders2 - derms^2)
    m <- anew[1]
    s <- anew[2]
    es <- exp(s)
    z1 <- ( y1 - m ) / es
    z2 <- ( y2 - m ) /es
    com <- sum(z1)
    ez1 <- exp(z1) 
    ez2 <- exp(z2)
    com1 <- sum(ez1)
    com2 <- sum(ez2)
  
    lik2 <- com1 - sum(ez1) - n1 * s - com2
    i <- 2
    while ( abs(lik2 - lik1) > tol ) {
      i <- i + 1
      lik1 <- lik2
    
      derm2 <-  - com1 / es^2 - com2 / es^2  
      derm <-  - n1/ es + - derm2 * es 
      ders <-  - com + sum(ez1 * z1) - n1 + sum(ez2 * z2)
      ders2 <-  com - sum(ez1 * z1^2) - sum(ez1 * z1) - sum(ez2 * z2^2) - sum(ez2 * z2)
      derms <- n1/es - sum(ez1 * z1) / es - com1/es - sum(ez2 * z2) / es - com2/es 
      anew <- c(m, s) - c(ders2 * derm - derms * ders, -derms * 
              derm + derm2 * ders)/(derm2 * ders2 - derms^2)
      m <- anew[1]
      s <- anew[2]
      es <- exp(s)
      z1 <- ( y1 - m ) / es
      z2 <- ( y2 - m ) /es
      com <- sum(z1)
      ez1 <- exp(z1) 
      ez2 <- exp(z2)
      com1 <- sum(ez1)
      com2 <- sum(ez2)
    
      lik2 <- com - com1 - n1 * s - com2
    }
    param <- c(exp(m), es)
    names(param) <- c("scale", "1/shape")
    mod <- list(iters = i, loglik = lik2 - sum(y1), param = param)
  } else  mod <- mod
  mod
} 


#[export]
gammapois.mle <- function(x, tol = 1e-07) {
  
  n <- length(x)
  #[export]
  slx <- sum( lgamma(x + 1) )
  sx <- sum(x) 
  m <- sx/n
  m2 <- sum(x^2)/n
  p <- 1 - m/(m2 - m^2)
  ini.ea <- max(0, m/p - m)
  eb <- ini.ea/m
  
  if (eb < 1) {
  lik1 <- sum( lgamma(x + ini.ea) ) - n * lgamma(ini.ea) + sx * log(eb/(1 + eb)) - n *ea * log1p(eb)

  dera <- sum( digamma(x + ini.ea) ) * ini.ea - n * digamma(ini.ea) * ea - n * ini.ea * log1p(eb)
  p <- eb / (1 + eb)
  derab <-  - n * ini.ea * p  
  derb <- sx * (1 - p) + derab
  dera2 <- dera + sum( trigamma(x + ini.ea) ) * ea^2 - n * trigamma(ini.ea) * ini.ea^2
  derb2 <-  - p * (1 - p) * sx - n * ini.ea * p * (1 - p)
  anew <- log(c(ini.ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
   
  ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
  lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)
  
  i <- 2 
  while ( lik2 - lik1 > tol ) {  
    i <- i + 1
    lik1 <- lik2
    dera <- sum( digamma(x + ea) ) * ea - n * digamma(ea) * ea - n * ea * log1p(eb)
    p <- eb / (1 + eb)
    derab <-  - n * ea * p  
    derb <- sx * (1 - p) + derab
    dera2 <- dera + sum( trigamma(x + ea) ) * ea^2 - n * trigamma(ea) * ea^2
    derb2 <-  - p * (1 - p) * sx - n * ea * p * (1 - p)
    anew <- log(c(ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
    ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
    lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)
  } 
  
  par <- c(ea, eb)
  names(par) <- c("shape", "scale")
  res <- list(iters = i, loglik = lik2, par = par) 

  } else  {
    ea <- ini.ea
    eb <- m/ea
    lik1 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n *ea * log1p(eb)

    dera <- sum( digamma(x + ea) ) * ea - n * digamma(ea) * ea - n * ea * log1p(eb)
    p <- eb / (1 + eb)
    derab <-  - n * ea * p  
    derb <- sx * (1 - p) + derab
    dera2 <- dera + sum( trigamma(x + ea) ) * ea^2 - n * trigamma(ea) * ea^2
    derb2 <-  - p * (1 - p) * sx - n * ea * p * (1 - p)
    anew <- log(c(ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
    ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
    lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)

    i <- 2 
    while ( lik2 - lik1 > tol ) {  
      i <- i + 1
      lik1 <- lik2
      dera <- sum( digamma(x + ea) ) * ea - n * digamma(ea) * ea - n * ea * log1p(eb)
      p <- eb / (1 + eb)
      derab <-  - n * ea * p  
      derb <- sx * (1 - p) + derab
      dera2 <- dera + sum( trigamma(x + ea) ) * ea^2 - n * trigamma(ea) * ea^2
      derb2 <-  - p * (1 - p) * sx - n * ea * p * (1 - p)
      anew <- log(c(ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
      ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
      lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)
    }
    par <- c(ea, eb)
    names(par) <- c("shape", "scale")
    res <- list(iters = i, loglik = lik2, par = par) 

  }  ##  end  if (eb < 1) {
  
  res
}

#gammapois.mle2 <- function(x, tol = 1e-07) {
  
#  n <- length(x)
#  slx <- sum( lgamma(x + 1) )
#  sx <- sum(x) 
#  m <- sx/n
#  m2 <- sum(x^2)/n
#  p <- 1 - m/(m2 - m^2)
#  ini.ea <- max(0, m/p - m)
#  eb <- ini.ea/m
  
#  lik1 <- sum( lgamma(x + ini.ea) ) - n * lgamma(ini.ea) + sx * log(eb/(1 + eb)) - n *ea * log1p(eb)

#  dera <- sum( digamma(x + ini.ea) ) * ini.ea - n * digamma(ini.ea) * ea - n * ini.ea * log1p(eb)
#  p <- eb / (1 + eb)
#  derab <-  - n * ini.ea * p  
#  derb <- sx * (1 - p) + derab
#  dera2 <- dera + sum( trigamma(x + ini.ea) ) * ea^2 - n * trigamma(ini.ea) * ini.ea^2
#  derb2 <-  - p * (1 - p) * sx - n * ini.ea * p * (1 - p)
#  anew <- log(c(ini.ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
   
#  ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
#  lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)
  
#  i <- 2 
#  while ( lik2 - lik1 > tol ) {  
#    i <- i + 1
#    lik1 <- lik2
#    dera <- sum( digamma(x + ea) ) * ea - n * digamma(ea) * ea - n * ea * log1p(eb)
#    p <- eb / (1 + eb)
#    derab <-  - n * ea * p  
#    derb <- sx * (1 - p) + derab
#    dera2 <- dera + sum( trigamma(x + ea) ) * ea^2 - n * trigamma(ea) * ea^2
#    derb2 <-  - p * (1 - p) * sx - n * ea * p * (1 - p)
#    anew <- log(c(ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
#    ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
#    lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)
#  } 
  
#  par <- c(ea, eb)
#  names(par) <- c("shape", "rate")
#  res <- list(iters = i, loglik = lik2, par = par) 

#  if (lik2 < lik1) {
#    ea <- ini.ea
#    eb <- m/ea
#    lik1 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n *ea * log1p(eb)

#    dera <- sum( digamma(x + ea) ) * ea - n * digamma(ea) * ea - n * ea * log1p(eb)
#    p <- eb / (1 + eb)
#    derab <-  - n * ea * p  
#    derb <- sx * (1 - p) + derab
#    dera2 <- dera + sum( trigamma(x + ea) ) * ea^2 - n * trigamma(ea) * ea^2
#    derb2 <-  - p * (1 - p) * sx - n * ea * p * (1 - p)
#    anew <- log(c(ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
#    ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
#    lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)

#    i <- 2 
#    while ( lik2 - lik1 > tol ) {  
#      i <- i + 1
#      lik1 <- lik2
#      dera <- sum( digamma(x + ea) ) * ea - n * digamma(ea) * ea - n * ea * log1p(eb)
#      p <- eb / (1 + eb)
#      derab <-  - n * ea * p  
#      derb <- sx * (1 - p) + derab
#      dera2 <- dera + sum( trigamma(x + ea) ) * ea^2 - n * trigamma(ea) * ea^2
#      derb2 <-  - p * (1 - p) * sx - n * ea * p * (1 - p)
#      anew <- log(c(ea, eb)) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
#      ea <- exp(anew[1])     ;   eb   <- exp(anew[2])
#      lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n * ea * log1p(eb)
#    }
#    par <- c(ea, eb)
#    names(par) <- c("shape", "scale")
#    res <- list(iters = i, loglik = lik2, par = par) 

#  }  ##  end  if (lik2 < lik1) {
  
#  res
#}
  
# fun = function(par, x, sx, n) {
#   ea <- exp(par[1])
#   eb <- exp(par[2])
#   lik2 <- sum( lgamma(x + ea) ) - n * lgamma(ea) + sx * log(eb/(1 + eb)) - n *ea * log1p(eb)
#   -lik2
# }
# optim( rnorm(2), fun, x = x, sx = sum(x), n= length(x) )


#[export]
halfcauchy.mle <- function(x, tol = 1e-07) {
   n <- length(x)
   es <- 0.5 * (Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4))
   logs <- log(es)
   x2 <- x^2
   down <- 1/(x2 + es^2) 
   lik1 <- n * logs + sum( log(down) )
   der <- n - 2 * es^2 * sum(down)
   der2 <-  - 4 * es^4 * sum(down^2)
   logs <- logs - der/der2
   es <- exp(logs)
   down <- 1/(x2 + es^2)   
   lik2 <- n * logs + sum( log(down) )
   i <- 2
   while ( lik2 - lik1 > tol ) {
     i <- i + 1
     lik1 <- lik2
     der <- n - 2 * es^2 * sum(down)
     der2 <-  - 4 * es^4 * sum(down^2)
     logs <- logs - der/der2
     es <- exp(logs)
     down <- 1/(x2 + es^2)   
     lik2 <- n * logs + sum( log(down) )
   } 
  list(iters = i, loglik = lik2 - n * log(2/pi), scale = es)
}


#[export]
cauchy0.mle<-function(x, tol = 1e-07) {
   n <- length(x)
   es <- 0.5 * (Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4))
   logs <- log(es)
   x2 <- x^2
   down <- 1/(x2 + es^2)   
   lik1 <- n * logs + sum( log(down) )
   der <- n - 2 * es^2 * sum(down)
   der2 <-  - 4 * es^4 * sum(down^2)
   logs <- logs - der/der2
   es <- exp(logs)
   down <- 1/(x2 + es^2)   
   lik2 <- n * logs + sum( log(down) )
   i <- 2
   while ( lik2 - lik1 > tol ) {
     i <- i + 1
     lik1 <- lik2
     der <- n - 2 * es^2 * sum(down)
     der2 <-  - 4 * es^4 * sum(down^2)
     logs <- logs - der/der2
     es <- exp(logs)
     down <- 1/(x2 + es^2)   
     lik2 <- n * logs + sum( log(down) )
   } 
  list(iters = i, loglik = lik2 - n * log(pi), scale = es)
}     


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

  
#[export]
powerlaw.mle <- function(x) {
  n <- length(x)
  x1 <- min(x)
  com <- sum( log(x) ) - n * log(x1)
  a <- 1 + n / com
  loglik <- n * log( (a - 1) / x1 ) - a * com 
  list(alpha = a, loglik = loglik)
}


#[export]
purka.mle <- function(x, tol = 1e-07) {
  if ( !is.matrix(x) )  x <- cbind( cos(x), sin(x) )
  p <- dim(x)[2]

  theta <- Rfast::mediandir(x)
  a <- x %*% theta
  a[ abs(a) > 1 ] <- 1
  A <- sum( acos(a) )
  n <- dim(x)[1]
  circle <- function(a, A, n)  n * log(a) - n * log(2) - n * log( 1 - exp( - a * pi ) ) - a * A
  sphere <- function(a, A, n)  n * log(a^2 + 1) - n * log(2 * pi) - n * log( 1 + exp( - a * pi ) ) - a * A
  hypersphere <- function(a, A, n) {
    n * lgamma(p/2) - 0.5 * n * p * log(pi) + n * ( log(besselI(a, p - 1, expon.scaled = TRUE)) + a ) - a * A
  }

  if (p == 2) {
    lika <- optimize(circle, c(0.001, 30000), maximum = TRUE, A = A, n = n, tol = tol) 
    a <- lika$maximum  ## estimated kappa
    f2 <-  -n / a^2 + n * pi^2 * exp( -a * pi)/ ( 1 - exp(-a * pi) )^2
  } else if (p == 3) {
    lika <- optimize(sphere, c(0.001, 30000), maximum = TRUE, A = A, n = n, tol = tol) 
    a <- lika$maximum  ## estimated kappa
    f2 <-  - (2 * a^2 * n - 2 * n) / (a^2 + 1)^2 - n * pi^2 * ( 1 + exp( a * pi) )^(-2) * exp( a * pi )    
  } else {
    lika <- optimize(hypersphere, c(0.001, 30000), maximum = TRUE, A = A, n = n, tol = tol)
    a <- lika$maximum  ## estimated kappa
    up1 <- ( besselI(a, p - 3) + 2 * besselI(p - 1, a) + besselI(p + 1, a) ) * besselI(p - 1, a)
    up2 <- ( besselI(a, p - 2) + 2 * besselI(p, a) )^2
    f2 <- 0.5 * n * ( up1 - up2 ) / besselI(p - 1, a)^2
  }
  ## f2 is the second derivative of the log-likelihood w.r.t alpha
    
  list( theta = theta, alpha = a, loglik = lika$objective, alpha.sd = 1 / sqrt( - f2 ) ) 
}


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

    
#[export]
trunccauchy.mle <- function (x, a, b, tol = 1e-07) {
    n <- length(x)
    m <- Rfast::med(x)
    es <- 0.5 * (Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4))
    logs <- log(es)
    y <- x - m
    y2 <- y^2
    lik1 <- n * logs - sum(log(es^2 + y2))
    down <- 1/(es^2 + y2)
    down2 <- down^2
    derm <- 2 * sum(y * down)
    ders <- n - 2 * es^2 * sum(down)
    derm2 <- 2 * sum((y2 - es^2) * down2)
    ders2 <- -2 * es^2 * (derm2 + 2 * es^2 * sum(down2))
    derms <- -4 * es^2 * sum(y * down2)
    m <- m - (ders2 * derm - derms * ders)/(derm2 * ders2 - derms^2)
    logs <- logs - (-derms * derm + derm2 * ders)/(derm2 * ders2 - 
        derms^2)
    y <- x - m
    y2 <- y^2
    es <- exp(logs)
    lik2 <- n * logs - sum(log(es^2 + y2))
    i <- 2
    while (lik2 - lik1 > tol) {
        i <- i + 1
        lik1 <- lik2
        down <- 1/(es^2 + y2)
        down2 <- down^2
        derm <- 2 * sum(y * down)
        ders <- n - 2 * es^2 * sum(down)
        derm2 <- 2 * sum((y2 - es^2) * down2)
        ders2 <- -2 * es^2 * (derm2 + 2 * es^2 * sum(down2))
        derms <- -4 * es^2 * sum(y * down2)
        m <- m - (ders2 * derm - derms * ders)/(derm2 * ders2 - 
            derms^2)
        logs <- logs - (-derms * derm + derm2 * ders)/(derm2 * 
            ders2 - derms^2)
        y <- x - m
        y2 <- y^2
        es <- exp(logs)
        lik2 <- n * logs - sum(log(es^2 + y2))
    }
    param <- c(m, es)
    names(param) <- c("location", "scale")
    tr <- atan(b) - atan(a)
    list(iters = i, loglik = lik2 - n * log(tr), param = param)
}


#[export]
truncexpmle <- function(x, b, tol = 1e-07) {
  trexp <- function(lam, sx, b, n) {
    - n * log(lam) - sx/lam - n * log( 1 - exp(-b/lam) )   
  }
  mod <- optimise(trexp, c(0, b), sx = sum(x), b = b, n = length(x), 
         tol = tol, maximum = TRUE )
  list(loglik = mod$objective, lambda = mod$minimum)
} 	


#[export]
zigamma.mle <- function(x, tol = 1e-07) {
  n <- length(x)
  x1 <- x[x > 0]
  n1 <- length(x1)
  n0 <- n - n1
  prob <- n1/n
  lik0 <- n0 * log(1 - prob) + n1 * log(prob)
  mod <- Rfast::gammamle(x1, tol = tol)
  param <- c(prob, mod$param)
  names(param) <- c("prop1", "shape", "scale")
  list(iters = mod$iters, loglik = sum(lik0, mod$loglik, na.rm = TRUE), param = param)
}


#[export]
zil.mle <- function(x) {

  n <- length(x)
  x1 <- x[x > 0]
  n1 <- length(x1)
  n0 <- n - n1
  prob <- n1/n
  lik0 <- n0 * log(1 - prob) + n1 * log(prob)

  lx1 <- log(x1)
  lx2 <- log(1 - x1)
  y <- lx1 - lx2
  sy <- sum(y)
  m <- sy/n1
  s <- (sum(y^2) - n1 * m^2)/n1
  loglik <- sum(dnorm(y, m, sqrt(s), log = TRUE)) - sy
  param <- c(prob, m, n1 * s/(n1 - 1))
  names(param) <- c("prop1", "mean", "unbiased variance")
  list(loglik = sum(lik0, loglik, na.rm = TRUE), param = param)
}


#[export]
ziweibull.mle <- function(x, tol = 1e-07) {
  n <- length(x)
  x1 <- x[x > 0]
  n1 <- length(x1)
  n0 <- n - n1
  prob <- n1/n
  lik0 <- n0 * log(1 - prob) + n1 * log(prob)
  mod <- Rfast::weibull.mle(x1, tol = tol)
  param <- c(prob, mod$param)
  names(param) <- c("prop1", "shape", "scale")
  list(iters = mod$iters, loglik = sum(lik0, mod$loglik, na.rm = TRUE), param = param)
}


#[export]
gnormal0.mle <- function(x, tol = 1e-06) {
  n <- length(x)
  xabs <- abs(x)

  fun <- function(b, xabs, n) {
    y <- xabs^b
    sy <- sum(y)
    1 + digamma(1/b) / b - sum( y * log(xabs) ) / sy + log( b/n * sy ) / b
  }
  
  b <- mean(xabs) / sqrt( mean(x^2) )
  mod <- uniroot(fun, lower = max(1e-5, b - 100 * b), upper = b + min(30, 200 * b), tol = tol, xabs = xabs, n = n )
  b <- mod$root
  a <- ( b/n * sum( xabs^b ) ) ^ ( 1/b )
  loglik <- n * log(b) - sum( xabs^b )/a^b - n * log(2 * a) - n * lgamma(1/b)        
  param <- c(a, b)
  names(param) <- c( "alpha", "beta" )
  list(iters = mod$iter, loglik = loglik, param = param)
}


#[export]
unitweibull.mle <- function(x, tol = 1e-07, maxiters = 100) {
  lx <-  - log(x)
  mod <- Rfast::weibull.mle( lx, tol = tol, maxiters = maxiters )
  param <- as.vector( mod$param )
  names(param) <- c("alpha", "beta")
  a <- mod$param[1]   ;   b <- mod$param[2]
  n <- length(x)
  loglik <- sum(lx) + n * log(a * b) + (b - 1) * sum( log(lx) ) - a * sum(lx^b)
  list(iters = mod$iters, loglik = loglik, param = param)
} 
