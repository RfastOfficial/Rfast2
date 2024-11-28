#[export]
lm.parboot <- function(x, y, R = 1000) {
   x <- model.matrix( y~., data.frame(x) )
   xx <- tcrossprod( Rfast::spdinv( crossprod(x) ), x ) 
   be <- xx %*% y
   dm <- dim(x)
   n <- dm[1]   ;    p <- dm[2]
   est <- as.vector( x %*% be ) 
   s <- Rfast::Var(y - est) * (n - 1)/(n - p)
   z <- Rfast::matrnorm(n, R)
   booty <- sqrt(s) * z + est
   xx %*% booty
}    


#[export]
lm.nonparboot <- function (x, y, R = 1000) {
    x <- model.matrix(y ~ ., data.frame(x))
    xx <- tcrossprod(Rfast::spdinv(crossprod(x)), x)
    be <- xx %*% y
    dm <- dim(x)
    n <- dm[1]
    p <- dm[2]
    est <- as.vector(x %*% be)
    res <- y - est
    s <- Rfast::Var(res) * (n - 1)/(n - p)
    z <- matrix(res, ncol = R)
    z <- Rfast::colShuffle(z)
    booty <- sqrt(s) * z + est
    xx %*% booty
}