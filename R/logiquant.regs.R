#[export]
logiquant.regs <- function(y, x, logged = FALSE) {
  m <- Rfast::med(y)
  y[y > m] <- 1
  y[y != 1] <- 0
  Rfast::univglms(y, x, oiko = "binomial", logged = logged)
}
  
