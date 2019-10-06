#[export]
colspml.mle <- function(x, tol = 1e-07, maxiters = 100, parallel = FALSE) {
   res <- .Call( Rfast2_colspml_mle,x, tol, maxiters, parallel)
   colnames(res) <- c("mu1", "mu2", "gamma", "loglik")
   res
}
