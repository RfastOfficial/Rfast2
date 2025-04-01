
#[export]
bessel <- function(x, nu, type = "I", expon.scaled = FALSE) {
  .Call(Rfast2_bessel, x, nu, type, expon.scaled)
}
