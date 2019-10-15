#[export]
halfcauchy.mle <- function(x, tol = 1e-07) {
   n <- length(x)
   es <- 0.5 * (Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4))
   logs <- log(es)
   x2 <- x^2
   down <- 1/(x2 + es^2)   
   lik1 <- n * logs + sum( log(down) )
   der <- n - 2 * es^2 * sum(down)
   der2 <-  - 4 * es^4 * sum(down^2)
   logs <- logs - der/der2
   es <- exp(logs)
   down <- 1/(x2 + es^2)   
   lik2 <- n * logs + sum( log(down) )
   i <- 2
   while ( lik2 - lik1 > tol ) {
     i <- i + 1
     lik1 <- lik2
     der <- n - 2 * es^2 * sum(down)
     der2 <-  - 4 * es^4 * sum(down^2)
     logs <- logs - der/der2
     es <- exp(logs)
     down <- 1/(x2 + es^2)   
     lik2 <- n * logs + sum( log(down) )
   } 
  list(iters = i, loglik = lik2 - n * log(2/pi), scale = es)
}


#[export]
cauchy0.mle <- function(x, tol = 1e-07) {
   n <- length(x)
   es <- 0.5 * (Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4))
   logs <- log(es)
   x2 <- x^2
   down <- 1/(x2 + es^2)   
   lik1 <- n * logs + sum( log(down) )
   der <- n - 2 * es^2 * sum(down)
   der2 <-  - 4 * es^4 * sum(down^2)
   logs <- logs - der/der2
   es <- exp(logs)
   down <- 1/(x2 + es^2)   
   lik2 <- n * logs + sum( log(down) )
   i <- 2
   while ( lik2 - lik1 > tol ) {
     i <- i + 1
     lik1 <- lik2
     der <- n - 2 * es^2 * sum(down)
     der2 <-  - 4 * es^4 * sum(down^2)
     logs <- logs - der/der2
     es <- exp(logs)
     down <- 1/(x2 + es^2)   
     lik2 <- n * logs + sum( log(down) )
   } 
  list(iters = i, loglik = lik2 - n * log(pi), scale = es)
}     


   