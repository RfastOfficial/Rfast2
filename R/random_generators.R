
#[export]
Runif <- function(n, min = 0, max = 1){
	.Call(Rfast2_Runif, n, min, max)
}

#[export]
Sample <- function(x, size, replace = FALSE){
	.Call(Rfast2_Sample, x, size, replace)
}

#[export]
Sample.int <- function(n, size = n, replace = FALSE){
	.Call(Rfast2_Sample_int, n, size, replace)
}

#[export]
Rbeta <- function(n, alpha, beta) {
    .Call(Rfast2_Rbeta, n, alpha, beta)
}

#[export]
Rexp <- function(n, rate = 1.0) {
    .Call(Rfast2_Rexp, n, rate)
}

#[export]
Rchisq <- function(n, df) {
    .Call(Rfast2_Rchisq, n, df)
}

#[export]
Rgamma <- function(n, shape, rate = 1.0) {
    .Call(Rfast2_Rgamma, n, shape, rate)
}

#[export]
Rgeom <- function(n, prob) {
    .Call(Rfast2_Rgeom, n, prob)
}

#[export]
Rcauchy <- function(n, location = 0, scale = 1) {
    .Call(Rfast2_Rcauchy, n, location, scale)
}

#[export]
Rt <- function(n, df, ncp = 0) {
    .Call(Rfast2_Rt, n, df, ncp)
}