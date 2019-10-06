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
