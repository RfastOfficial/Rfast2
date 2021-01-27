#[export]
fp <- function(y, x, a, di = NULL, type = "logistic", full = FALSE, 
               seb = FALSE, tol = 1e-07, maxiters = 100) {
			   
			   
  if (type == "normal") {
    fun <- function(a, y, x)  abs( cor(y, x^a) ) )
    mod <- optimize(fun, a, y = y, x = x, maximum = TRUE)
    a <- mod$maximum
    mod <- Rfast::lmfit(cbind(1, x^a), y)  
   
  } else if (type == "logistic")
    fun <- function(a, y, x)  Rfast::glm_logistic(y, x^a, tol = tol, maxiters = maxiters)$devi
    mod <- optimize(fun, a, y = y, x = x)
    a <- mod$minimum
    mod <- Rfast::glm_logistic(y, x^a, full = full, tol = tol, maxiters = maxiters)

  } else if ( type == "poisson" )
    fun <- function(a, y, x)  Rfast::glm_poisson(y, x^a, tol = tol)$devi
    mod <- optimize(fun, a, y = y, x = x)
    a <- mod$minimum
    mod <- Rfast::glm_poisson(y, x^a, full = full, tol = tol)

  } else if ( type == "spml" )
    fun <- function(a, y, x)  Rfast::spml.reg(y, x^a, tol = tol, maxiters = maxiters)$loglik
    mod <- optimize(fun, a, y = y, x = x, maximum = TRUE)
    a <- mod$maximum
    mod <- Rfast::spml.reg(y, x^a, tol = tol, seb =seb, maxiters = maxiters)

  } else if ( type == "gamma" )
    fun <- function(a, y, x)  Rfast2::gammareg(y, x^a, tol = tol, maxiters = maxiters)$deviance
    mod <- optimize(fun, a, y = y, x = x)
    a <- mod$minimum
    mod <- Rfast2::gammareg(y, x^a, tol = tol, maxiters = maxiters)
  
  } else if ( type == "normlog" )
    fun <- function(a, y, x)  Rfast::normlog.reg(y, x^a, tol = tol, maxiters = maxiters)$deviance
    mod <- optimize(fun, a, y = y, x = x)
    a <- mod$minimum
    mod <- Rfast::normlog.reg(y, x^a, tol = tol, maxiters = maxiters)

  } else if ( type == "weibull" )
    fun <- function(a, y, x)  Rfast2::censweib.reg(y, x^a, di = di, tol = tol, maxiters = maxiters)$loglik
    mod <- optimize(fun, a, y = y, x = x, maximum = TRUE)
    a <- mod$maximum
    mod <- Rfast2::censweib.reg(y, x^a, di, tol = tol, maxiters = maxiters) 

  } else if ( type == "negbin" )
    fun <- function(a, y, x)  Rfast2::negbin.reg(y, x^a, tol = tol, maxiters = maxiters)$loglik
    mod <- optimize(fun, a, y = y, x = x, maximum = TRUE)
    a <- mod$maximum
    mod <- Rfast2::negbin.reg(y, x^a, tol = tol, maxiters = maxiters) 
  }

  list(a = a, mod = mod)
}













