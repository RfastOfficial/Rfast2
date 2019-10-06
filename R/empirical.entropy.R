#[export]
empirical.entropy <- function(x, k = NULL, pretty = FALSE) {
  if ( is.null(k) ) {
    n <- length(x)
    h <- 2 * diff( Rfast2::Quantile( x, probs = c(0.25, 0.75) ) ) / n^(1/3)
    mm <- Rfast::min_max(x)
    k <- ceiling( ( mm[2] - mm[1] )/ h ) 
  }
  if ( pretty )  {
    a <- pretty(x, n = k + 1) 
  } else  a <- seq( mm[1] - 1, mm[2] + 1, length = k + 1 )
  freqs <- Rfast::Table( cut(x, a) )
   - sum( freqs * log(freqs), na.rm = TRUE )
}