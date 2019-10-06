#[export]
welch.tests <- function(y, x, logged = FALSE, parallel = FALSE) {
   res <- .Call( Rfast2_welch_tests,x, y, logged, parallel)
   colnames(res) <- c("F stat", "pvalue")
   res
}




