#[export]
big.knn <- function(xnew, y, x, k = 2:100, type = "R") {

  if ( !is.matrix(xnew) )  xnew <- matrix(xnew, nrow = 1)
  di <- RANN::nn2( data = x, query = xnew, k = max(k) )$nn.idx
  nu <- dim(xnew)[1]
  denom <- 1:p
  nk <- length(k)
  est <- matrix(nrow = nu, ncol = nk + 1)
  for ( i in 1:nu ) est[i, ] <- cumsum( y[ di[i, ] ] ) / denom
  est <- est[, -c( min(k) - 1 ) ]
  colnames(est) <- paste("k=", k, sep = "")
  est
}
