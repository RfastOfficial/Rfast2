#[export]
jack.mean <- function(x) {
  n <- length(x)
  sum(x) * (n - 1) / n^2 
}



#[export]
coljack.means <- function(x) {
  n <- dim(x)[1]
  Rfast::colsums(x) * (n - 1) / n^2
}



#[export]
rowjack.means <- function(x) {
  n <- dim(x)[2]
  Rfast::rowsums(x) * (n - 1) / n^2
}



