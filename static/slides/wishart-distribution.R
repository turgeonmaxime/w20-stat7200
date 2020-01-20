## ----setup, include=FALSE---------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ---------------------------------------------------------------
B <- 1000
n <- 10; p <- 4

traces <- replicate(B, {
  Z <- matrix(rnorm(n*p), ncol = p)
  W <- crossprod(Z)
  sum(diag(W))
})


## ---------------------------------------------------------------
hist(traces, 50, freq = FALSE)
lines(density(rchisq(B, df = n*p)))

