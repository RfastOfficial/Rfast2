colbeta.mle <- function (x, tol = 1e-07, maxiters = 100, parallel = FALSE {
    res <- .Call(Rfast2_colbeta_mle, x, tol, parallel, maxiters)
    colnames(res) <- c("alpha", "betag", "loglik")
    res
}