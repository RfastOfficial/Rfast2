#[export]
multinom.reg <- function(y, x, tol = 1e-07, maxiters = 100) {
  x <- model.matrix( y ~. , data.frame(x) )
  .Call( Rfast2_multinom_reg,y, x, tol, maxiters)
}  




