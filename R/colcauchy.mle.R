#[export]
colcauchy.mle <- function (x, tol = 1e-07, maxiters = 100, parallel = FALSE) {
    res <- .Call(Rfast2_colcauchy_mle, x, tol, parallel, maxiters)
    colnames(res) <- c("loglik", "location", "scale")
    res
}
