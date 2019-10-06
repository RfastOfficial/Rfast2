#[export]
gumbel.reg <- function(y, x, tol = 1e-07, maxiters = 100) {
   
  X <- model.matrix(y ~., data.frame(x) )
  sx <- Rfast::colsums(X)
  dm <- dim(X)
  n <- dm[1]
  p <- dm[2] 
  mod <- Rfast::lmfit(X, y)
  be <- mod$be
  s <- sqrt( sum( mod$residuals^2 )/(n - p) )
  z <- as.vector( mod$residuals ) / s 
  exp_z <- exp(-z)  
  lik1 <-  - n * log(s) - sum(z) - sum(exp_z)
  com <- exp_z * X
  ders <-  - n + sum(z) - sum(exp_z * z)
  ders2 <-  - n - ders - sum(exp_z * z^2)
  derb <- sx/s - Rfast::colsums(com) / s
  derb2 <-  - crossprod(com, X) / s^2
  ## derbs <-  - sx/s - Rfast::colsums(com * z) / s^2 + Rfast::colsums(com)/s^2 
  be <- be - solve(derb2, derb)
  s <- log(s) - ders/ders2 
  s <- exp(s)
  z <- as.vector(y - X %*% be) / s
  exp_z <- exp(-z) 
  lik2 <-  - n * log(s) - sum(z) - sum(exp_z)
  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters) {
    i <- i + 1
    lik1 <- lik2
    com <- exp_z * X
    ders <-  - n + sum(z) - sum(exp_z * z)
    ders2 <-  - n - ders - sum(exp_z * z^2)
    derb <- sx/s - Rfast::colsums(com)/s
    derb2 <-  - crossprod(com, X) / s^2
    ## derbs <-  - sx/s - Rfast::colsums(com * z)/s^2 + Rfast::colsums(com)/s^2
    be <- be - solve(derb2, derb)
    s <- log(s) - ders/ders2 
    s <- exp(s)
    z <- as.vector(y - X %*% be) / s
    exp_z <- exp(-z)  
    lik2 <-  - n * log(s) - sum(z) - sum(exp_z)
  }
  list(be = be, sigma = s, loglik = lik2, iters = i)
}




#gumb <- function(y, x) {
#  X <- model.matrix(y~., data.frame(x) )
#  mod <- lm.fit(X, y)
#  dm <- dim(X)[2] 
#  be <- coef(mod)
#  s <- sqrt( sum( mod$residuals^2 )/mod[[ 8 ]] )
#  pa <- c(coef(mod), log(s) ) 
#  n <- length(y)   
#  gumbel <- function(pa) {
#    be <- pa[1:dm]   ;  s <- exp( pa[dm + 1] )
#    m <- X %*% be
#    z <- y - m
#    n * log(s) + sum(z)/s + sum( exp(-z/s) )  
#  }
#  nlm(gumbel, pa)
#}



