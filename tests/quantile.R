x <- c(9, 1, 5, 3)
probs <- c(1, 0, 0.75, 0.25, 0.5)
expected <- unname(stats::quantile(x, probs, type = 7))
stopifnot(isTRUE(all.equal(Rfast2::Quantile(x, probs),
                           expected)))

x <- 16:1
probs <- c(1, 0.35, 0)
expected <- unname(stats::quantile(x, probs, type = 7))
stopifnot(isTRUE(all.equal(Rfast2::Quantile(x, probs),
                           expected)))
