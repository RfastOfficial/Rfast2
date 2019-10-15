#[export]
depth.mahala <- function(x, data) {
  a <- 1 + Rfast::mahala(x, Rfast::colmeans(data), Rfast::cova(data) )
  1/a
}