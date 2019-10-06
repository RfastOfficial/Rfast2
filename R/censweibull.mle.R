#[export]
censweibull.mle <- function(x, di, tol = 1e-07) {
  y <- log(x)
  y1 <- y[di == 1]
  y2 <- y[di == 0]
  n <- length(x)
  n1 <- length(y1)
  n2 <- n - n1
  mod <- weibull.mle(x[di==1])
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
  
