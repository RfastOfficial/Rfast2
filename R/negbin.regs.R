negbin.regs <- function (y, x, tol = 1e-07, maxiters = 100, parallel = FALSE)
{
  mod <- .Call(Rfast2_negbin_regs, y, x, tol, maxiters, parallel)
  names(mod$info) <- c("iters", "BIC", "log-likelihood",
                       "dispersion")
  row.names(mod$info) <- names(x)

  names(mod$be) <- c("(intercept)","be")
  row.names(mod$be) <- names(x)
  list(info = mod$info, be = mod$be)
}
