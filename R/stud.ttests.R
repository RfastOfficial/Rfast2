#[export]
stud.ttests <- function(x, y = NULL, ina, logged = FALSE, parallel = FALSE) {

  if ( is.null(y) ) {
    x1 <- x[ ina == 1, ]
    x2 <- x[ ina == 2, ]
    n1 <- sum( ina == 1 )
    n2 <- length(ina) - n1
  } else {
    x1 <- x     ;    n1 <- dim(x1)[1]
    x2 <- y     ;    n2 <- dim(x2)[1]
  }

  m1 <- Rfast::colmeans(x1, parallel = parallel)
  m2 <- Rfast::colmeans(x2, parallel = parallel)
  s1 <- Rfast::colVars(x1, suma = n1 * m1, parallel = parallel) 
  s2 <- Rfast::colVars(x2, suma = n2 * m2, parallel = parallel)
  sp <- ( (n1 - 1) * s1 + (n2 - 1) * s2 ) / (n1 + n2 - 2) 
  stat <- ( m1 - m2 ) / sqrt(sp * (1/n1 + 1/n2) )
  dof <- n1 + n2 - 2
  if ( logged ) {
    pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE )
  } else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE )  
  result <- cbind(stat, pvalue, dof)

  result
}
