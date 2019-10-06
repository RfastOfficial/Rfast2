#[export]
fbed.reg <- function (y, x, alpha = 0.05, type = "logistic", K = 0, 
    backward = FALSE, parallel = FALSE, tol = 1e-07, maxiters = 100) {
    mod <- .Call(Rfast2_fbed_reg, y, x, alpha, type, id = integer(0), 
        K, backward, tol, parallel, maxiters)
    mod$ini <- as.vector(mod$startmod)
    mod$startmod <- NULL
    mod$res <- mod$colsfound
    mod$colsfound <- NULL
    if (dim(mod$res)[2] == 1) 
        mod$res <- matrix(0, 1, 2)
    colnames(mod$res) <- c("Vars", "stat")
    mod$info <- mod$kmatrix[, -1, drop = FALSE]
    rownames(mod$info) <- paste("K=", 0:K, sep = "")
    colnames(mod$info) <- c("Number of vars", "Number of tests")
    mod
}
