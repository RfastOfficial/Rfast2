#[export]
moranI <- function(x, w, scaled = FALSE, R = 999) {
  if ( !scaled )  w <- w / Rfast::rowsums(w)
  y <- w %*% x
  if (R > 1) {
    b <- Rfast::permcor(y, x, R = R) 
    res <- c( b[1] * Rfast::Var(y, std = TRUE) / Rfast::Var(x, std = TRUE), b[2] )
  } else  res <- c( cor(y, x,) * Rfast::Var(y, std = TRUE) / Rfast::Var(x, std = TRUE), NA )	
  names(res) <- c("Moran's I", "permutation p-value")
  res
}
