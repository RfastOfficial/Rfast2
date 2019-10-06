#[export]
km <- function(ti, di) {
  x <- cbind(ti, di)
  x <- x[order(ti), ]
  n <- dim(x)[1]
  id <- 1:n
  di <- x[, 2]
  ri <- n + 1 - id[ di > 0 ]
  time <- x[di>0, 1]
  event <- Rfast::Table(time, names = FALSE)
  mi <- cbind( time, ri, event )
  survival <- 1 - event/ri  
  survival <- cumprod(survival)
  cbind(mi, survival)
}

  
