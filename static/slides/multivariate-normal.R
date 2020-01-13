## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## ---- message = FALSE---------------------------------------------------------
set.seed(123)

n <- 1000; p <- 2
Z <- matrix(rnorm(n*p), ncol = p)

mu <- c(1, 2)
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
L <- t(chol(Sigma))


## ---- message = FALSE---------------------------------------------------------
Y <- L %*% t(Z) + mu
Y <- t(Y)

colMeans(Y)
cov(Y)

library(tidyverse)
Y %>% 
  data.frame() %>% 
  ggplot(aes(X1, X2)) +
  geom_density_2d()


## -----------------------------------------------------------------------------
library(mvtnorm)

Y <- rmvnorm(n, mean = mu, sigma = Sigma)

colMeans(Y)
cov(Y)

Y %>% 
  data.frame() %>% 
  ggplot(aes(X1, X2)) +
  geom_density_2d()

