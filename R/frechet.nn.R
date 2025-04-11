#[export]
frechet.nn <- function(x, di, a, k, parallel = FALSE, cores = 0) {
	.Call(Rfast2_frechet_nn, x, di, a, k, parallel, cores)
}