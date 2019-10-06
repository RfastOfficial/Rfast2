#[export]
covar <- function(y, x) {
 denom <- dim(x)[1] - 1
 ( Rfast::eachcol.apply(x, y) - Rfast::colmeans(x) * sum(y) ) / denom
}