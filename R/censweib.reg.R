
#[export]
censweib.reg <- function (y, x, di = NULL, tol = 1e-07, maxiters = 100) {
    X <- model.matrix(y ~ ., data.frame(x))
	mod<-NULL
	if ( is.null(di) ) {
	  mod <- .Call(Rfast2_censweib_reg, y, X, tol = tol, maxiters =  maxiters)
	} else {
	  mod <- .Call(Rfast2_censweib_reg, y, X, di, tol, maxiters)
	}
    rownames(mod$be) <- colnames(X)
    list(iters = mod$iters, loglik = mod$loglik, shape = mod$shape, 
        be = mod$be)
}
