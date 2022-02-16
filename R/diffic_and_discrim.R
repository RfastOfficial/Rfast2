#[export]
diffic <- function(x) {
  Rfast::colmeans(x)
}



#[export]
discrim <- function(x, frac = 1/3) {
  n <- dim(x)[1]
  N <- floor(n * frac)
  score <- Rfast::rowmeans(x)  
  score <- x[order(score), ]
  up <- x[(n + 1 - N):n, ]
  low <- x[1:N, ]
  up.score <- Rfast::colsums(up)
  low.score <- Rfast::colsums(low)
  a <- (up.score - low.score)/N
  names(a) <- colnames(x)
  a
} 