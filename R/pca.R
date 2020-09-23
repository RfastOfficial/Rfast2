#[export]
pca <- function(x, center = TRUE, scale = TRUE, k = NULL, vectors = FALSE) {
   dm <- dim(x)
   n <- dm[1]
   p <- dm[2]
   if ( center | scale )   x <- Rfast::standardise(x, center = center, scale = scale)
   if (n < p ) {
     res <- Rfast::hd.eigen(x, center = FALSE, scale = FALSE, k = k, vectors = vectors)
   } else {
     if ( is.null(k) )  k <- p
     if ( !vectors ) {
        values <- svd(x, nu = 0, nv = 0)$d^2 / (n - 1)
        vectors <- NULL
     } else {
       xsvd <- svd( x, nu = 0)
       values <- xsvd$d^2 / (n - 1)
       vectors <- xsvd$v[, 1:k, drop = FALSE]
     }
     res <- list(values = values[1:k], vectors = vectors)
   }
   res
}
     
   