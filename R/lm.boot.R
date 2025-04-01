#[export]
lm.boot <- function(x, y, R = 1000, type = "param") {
   x <- model.matrix( y~., data.frame(x) )
   xx <- tcrossprod( Rfast::spdinv( crossprod(x) ), x ) 
   be <- xx %*% y
   dm <- dim(x)
   n <- dm[1]   ;    p <- dm[2]
   est <- as.vector( x %*% be ) 
   res <- y - est
   s <- Rfast::Var(res) * (n - 1) / (n - p)
    
   if ( type == "param" ) {
     z <- Rfast::matrnorm(n, R)
   } else if ( type == "nonparam" ) {
     z <- Rfast2::Sample(res, n * R, replace = TRUE)
     dim(z) <- c(n, R) 
   }	 
   booty <- sqrt(s) * z + est
   xx %*% booty
}    

#[export]
lm.parboot <- function(x, y, R = 1000) {
   .Deprecated(new = "lm.boot(type = \"param\")", old = "lm.parboot")
}    


#[export]
lm.nonparboot <- function (x, y, R = 1000) {
   .Deprecated(new = "lm.boot(type = \"nonparam\")", old = "lm.nonparboot")
}