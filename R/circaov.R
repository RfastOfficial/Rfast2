#[export]
hcf.circaov <- function (u, ina)  {
  n <- length(u)
  ina <- as.numeric(ina)
  g <- max(ina)
  x1 <- cos(u)
  x2 <- sin(u)
  Ci <- Rfast::group(x1, ina)
  Si <- Rfast::group(x2, ina)
  Ri <- sqrt(Ci^2 + Si^2)
  V <- sum(Ri)
  C <- sum(Ci)
  S <- sum(Si)
  R <- sqrt(C^2 + S^2)

  mu <- atan(S/C) + pi * ( C < 0 )
  con <- sum(cos(u - mu))
  k1 <- (1.28 - 0.53 * R^2/n^2) * tan(0.5 * pi * R/n)
  if (k1 < 710) {
    der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
    a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
    der2 <- n * a/besselI(k1, 0)^2
    k2 <- k1 + der/der2
    while (abs(k1 - k2) > 1e-8) {
      k1 <- k2
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
    }
  } else k2 <- k1
  kapa <- k2

  if (kapa > 2) {
    Ft <- (n - g) * (V - R)/(g - 1)/(n - V)
    pvalue <- pf(Ft, g - 1, n - g, lower.tail = FALSE)
  } else if (kapa < 2 & kapa > 1) {
    Ft <- (1 + 3/(8 * kapa)) * (n - g) * (V - R) / ( (g - 1) * (n - V) )
    pvalue <- pf(Ft, g - 1, n - g, lower.tail = FALSE)
  } else {
    Ft <- NA
    pvalue <- NA
  }
  res <- c(Ft, pvalue, kapa)
  names(res) <- c("test", "p-value", "kappa")
  res
}


#[export]
het.circaov <- function(u, ina) {
  n <- length(u)
  ina <- as.numeric(ina)
  g <- max(ina)
  ni <- tabulate(ina)
  kapa <- numeric(g)
  x1 <- cos(u)
  x2 <- sin(u)
  C <- Rfast::group(x1, ina)
  S <- Rfast::group(x2, ina)
  mi <- atan(S/C) + pi * as.numeric(C < 0)
  Ri <- sqrt(C^2 + S^2)

  ki <- (1.28 - 0.53 * Ri^2/ni^2) * tan(0.5 * pi * Ri/ni)
  coni <- ki
  for (i in 1:g) {
    n <- ni[i]
    coni[i] <- sum( cos( u[ina == i] - mi[i] ) )
    con <- coni[i]
    k1 <- ki[i]
    if (k1 < 710) {
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
      while ( abs(k2 - k1) > 1e-08 ) {
        k1 <- k2
        der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
      }
    } else k2 <- k1
    kapa[i] <- k2
  }

  Rw <- sqrt(sum(kapa * Ri * cos(mi))^2 + sum( kapa * Ri * sin(mi))^2 )
  Ta <- 2 * (sum(kapa * Ri) - Rw)
  pvalue <- pchisq(Ta, g - 1, lower.tail = FALSE)
  res <- c(Ta, pvalue)
  names(res) <- c("test", "p-value")
  res
}


#[export]
lr.circaov <- function(u, ina) {
  ina <- as.numeric(ina)
  g <- max(ina)
  x <- cbind(cos(u), sin(u))
  rsi <- rowsum(x, ina)
  Ri <- sqrt(Rfast::rowsums(rsi^2))
  ni <- tabulate(ina)
  mi <- rsi/ni
  mi <- mi/sqrt(Rfast::rowsums(mi^2))
  m <- Rfast::colmeans(x)
  m <- m/sqrt(sum(m^2))
  m <- matrix(rep(m, g), nrow = g, byrow = TRUE)

  n <- dim(x)[1]
  rs <- Rfast::colsums(rsi)
  mu <- atan(rs[2]/rs[1]) + pi * (rs[1] < 0)
  con <- sum( cos(u - mu) )
  R <- sqrt( sum(rs^2) )
  k1 <- (1.28 - 0.53 * R^2/n^2) * tan(0.5 * pi * R/n)
  if (k1 < 710) {
        der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1,
            0, expon.scaled = TRUE)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1,
            0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
        while (abs(k1 - k2) > 1e-08) {
            k1 <- k2
            der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1,
                0, expon.scaled = TRUE)
            a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1,
                0)/2 - besselI(k1, 1)^2
            der2 <- n * a/besselI(k1, 0)^2
            k2 <- k1 + der/der2
        }
    } else k2 <- k1
  kapa <- k2

  w <- kapa * sum( Ri * Rfast::rowsums( (mi - m)^2 ) )
  pvalue <- pchisq(w, g - 1, lower.tail = FALSE)
  res <- c(w, pvalue, kapa)
  names(res) <- c("test", "p-value", "kappa")
  res
}


#[export]
embed.circaov <- function (u, ina) {
  n <- length(u)
  ina <- as.numeric(ina)
  ni <- tabulate(ina)
  g <- max(ina)
  x1 <- cos(u)
  x2 <- sin(u)
  Ci <- Rfast::group(x1, ina)
  Si <- Rfast::group(x2, ina)
  Rbi <- (Ci^2 + Si^2)/ni^2
  C <- sum(Ci)
  S <- sum(Si)
  Rbar <- sqrt(C^2 + S^2)/n

  mu <- atan(S/C) + pi * (C < 0)
  con <- sum( cos(u - mu) )
  k1 <- (1.28 - 0.53 * Rbar^2) * tan(0.5 * pi * Rbar)
  if (k1 < 710) {
    der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
    a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
    der2 <- n * a/besselI(k1, 0)^2
    k2 <- k1 + der/der2
    while (abs(k1 - k2) > 1e-08) {
      k1 <- k2
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
    }
  } else k2 <- k1
  kapa <- k2

  Fb <- ( (sum(ni * Rbi) - n * Rbar^2)/(g - 1) ) / ( (n - sum(ni * Rbi) ) / (n - g) )
  Fc <- (1 - 1/(5 * kapa) - 1/(10 * kapa^2)) * Fb
  pvalue <- pf(Fc, g - 1, n - g, lower.tail = FALSE)
  res <- c(Fc, pvalue, kapa)
  names(res) <- c("test", "p-value", "kapa")
  res
}

