#[export]
negbin.reg <- function(y, x, tol = 1e-07, maxiters = 100) {
  x <- model.matrix( y ~. , data.frame(x) )
  mod <- .Call( Rfast2_negbin_reg,y, x, tol, maxiters)
  names(mod$info) <- c( "iters", "BIC", "log-likelihood", "dispersion" )
  list(info = mod$info, be = mod$be)
}  




