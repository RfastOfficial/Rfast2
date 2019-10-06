#[export]
weib.regs <- function (y, x, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100) {
   mod <- .Call( Rfast2_weib_regs,y, x, tol, logged, maxiters, parallel)
    colnames(mod) <- c("stat", "pvalue")
    mod
}
