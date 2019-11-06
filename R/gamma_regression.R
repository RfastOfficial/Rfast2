#[export]
gammareg <- function(y, x, tol = 1e-07, maxiters = 100) {
  mod <- Rfast::gammacon(y)
  gamma_reg(Y = y, X = x, mod = mod, tol = tol, maxiters = maxiters) 
}


#[export]
gammaregs <- function(y, x, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100) {
  gamma_regs(Y = y, X = x, tol = tol, logged = logged, parallel = (parallel > 0), maxiters = maxiters)
}
