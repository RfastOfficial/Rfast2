#[export]
leverage <- function(x) {
  Rfast::mahala(x, numeric(dim(x)[2]), crossprod(x) ) 
} 
