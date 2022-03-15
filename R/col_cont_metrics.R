#[export]
colmses <- function(y, yhat, parallel = FALSE) {
  Rfast::colmeans( (y - yhat)^2, parallel = parallel )
}


#[export]
colmaes <- function(y, yhat, parallel = FALSE) {
  Rfast::colmeans( abs(y - yhat), parallel = parallel )
}


#[export]
colpkl <- function(y, yhat, parallel = FALSE) {
 y1 <- 1 - y
 Rfast::colsums( y * log(y / yhat), na.rm = TRUE, parallel = parallel ) +  
 Rfast::colsums( y1 * log( y1 / (1 - yhat) ), na.rm = TRUE, parallel = parallel )
}


#[export]
colukl <- function(y, yhat, parallel = FALSE) {
 Rfast::colsums( y * log(y / yhat), na.rm = TRUE, parallel = parallel )
}

