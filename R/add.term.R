#[export]
add.term <- function (y, xinc, xout, devi_0, type = "logistic", logged = FALSE, 
    tol = 1e-07, maxiters = 100, parallel = FALSE) {
    .Call(Rfast2_add_term, y, cbind(1, xinc), xout, devi_0, type, tol, 
        logged, parallel, maxiters)
}
