#[export]
hp.reg <- function(y, x, full = FALSE, tol = 1e-07, maxiters = 100) {
  y1 <- y
  id <- which(y > 0)
  y1[id] <- 1
  prob <- Rfast::glm_logistic(x, y1, full = full, tol = tol, maxiters = maxiters)
  x <- model.matrix(y ~ ., data.frame(x) )
  mod <- Rfast2::ztp.reg(y[id], x[id, ], full = full, tol = tol, maxiters = maxiters)
  list(prob = prob, mod = mod)
}
  

    