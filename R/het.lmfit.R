#[export]
het.lmfit <- function(x, y, type = 1) {
   if ( type == 1 )  {
     xxinvtx <- solve( crossprod(x), t(x) )
     be <- xxinvtx %*% y
     u <- y - x %*% be
     be <- xxinvtx %*% log(u^2)  
     gi <- as.vector( x %*% be )
   } else if ( type == 2 ) {
     yhat <- x %*% be
     u <- y - yhat
     z <- cbind(1, yhat, yhat^2)  
     be <- solve( crossprod(z), crossprod(z, log(u^2)) )
     gi <- as.vector( z %*% be )
   }
   hi <- 1/sqrt( exp(gi) )
   be <- solve(crossprod(x, hi * x), crossprod(x, hi * y))
   e <- y - x %*% be
   list(be = be, residuals = e)
}
