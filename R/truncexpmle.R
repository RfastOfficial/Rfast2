#[export]
truncexpmle <- function(x, b, tol = 1e-07) {
  
  trexp <- function(lam, sx, b, n) {
    - n * log(lam) - sx/lam - n * log( 1 - exp(-b/lam) )   
  }
  mod <- optimise(trexp, c(0, b), sx = sum(x), b = b, n = length(x), 
         tol = tol, maximum = TRUE )
  list(loglik = mod$objective, lambda = mod$minimum)
} 