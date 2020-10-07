#[export]
censweib.reg <- function (y, x, di, tol = 1e-07, maxiters = 100) {
    X <- model.matrix(y ~ ., data.frame(x))
	if ( is.null(di) ) {
	  mod <- Rfast::weib.reg(y, X, tol = tol, maxiters =  maxiters)
	} else {
	  mod <- .Call(Rfast2_censweib_reg, y, X, di, tol, maxiters)
	}
    names(mod$be) <- colnames(X)
    list(iters = mod$iters, loglik = mod$loglik, shape = mod$shape, 
        be = mod$be)
}
