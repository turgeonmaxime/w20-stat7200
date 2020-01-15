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


## -----------------------------------------------------------------------------
# Ramus data, Timm (2002)
main_page <- "https://maxturgeon.ca/w20-stat7200/"
ramus <- read.csv(paste0(main_page, "Ramus.csv"))
head(ramus, n = 5)


## -----------------------------------------------------------------------------
var_names <- c("Age8", "Age8.5",
               "Age9", "Age9.5")

par(mfrow = c(2, 2))
for (var in var_names) {
  qqnorm(ramus[, var], main = var)
  qqline(ramus[, var])
}


## -----------------------------------------------------------------------------
ramus <- ramus[,var_names]
sigma_hat <- cov(ramus)

ramus_cent <- scale(ramus, center = TRUE, 
                    scale = FALSE)

D_vect <- apply(ramus_cent, 1, function(row) {
  t(row) %*% solve(sigma_hat) %*% row
})


## -----------------------------------------------------------------------------
qqplot(qchisq(ppoints(D_vect), df = 4),
       D_vect, xlab = "Theoretical Quantiles")
qqline(D_vect, distribution = function(p) {
  qchisq(p, df = 4)
  })


## -----------------------------------------------------------------------------
library(mvtnorm)
set.seed(123)

n <- 50; p <- 2

mu <- c(1, 2)
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = p)

Y <- rmvnorm(n, mean = mu, sigma = Sigma)


## -----------------------------------------------------------------------------
loglik <- function(mu, sigma, data = Y) {
  # Compute quantities
  y_bar <- colMeans(Y)
  quad_form <- t(y_bar - mu) %*% solve(sigma) %*%
    (y_bar - mu)
  
  -0.5*n*log(det(sigma)) -
    0.5*(n - 1)*sum(diag(solve(sigma) %*% cov(Y))) -
    0.5*n*drop(quad_form)
}


## -----------------------------------------------------------------------------
grid_xy <- expand.grid(seq(0, 2, length.out = 32), 
                       seq(0, 4, length.out = 32))

head(grid_xy, n = 5)


## -----------------------------------------------------------------------------
contours <- purrr::map_df(seq_len(nrow(grid_xy)), 
                          function(i) {
  # Where we will evaluate loglik
  mu_obs <- as.numeric(grid_xy[i,])
  # Evaluate at the pop covariance
  z <- loglik(mu_obs, sigma = Sigma)
  # Output data.frame
  data.frame(x = mu_obs[1],
             y = mu_obs[2],
             z = z)
})


## ---- message = FALSE---------------------------------------------------------
library(tidyverse)
library(ggrepel)
# Create df with pop and sample means
data_means <- data.frame(x = c(mu[1], mean(Y[,1])),
                         y = c(mu[2], mean(Y[,2])),
                         label = c("Pop.", "Sample"))


## ---- message = FALSE---------------------------------------------------------
ggplot(contours, aes(x, y)) + 
  geom_contour(aes(z = z)) + 
  geom_point(data = data_means) +
  geom_label_repel(data = data_means,
                   aes(label = label))


## -----------------------------------------------------------------------------
library(scatterplot3d)
with(contours, scatterplot3d(x, y, z))


## ---- echo = FALSE, eval = FALSE----------------------------------------------
## x <- seq(0, 2, length.out = 32)
## y <- seq(0, 4, length.out = 32)
## z <- matrix(NA, ncol = 32, nrow = 32)
## for (i in seq_len(32)) {
##   for (j in seq_len(32))
##     z[i,j] <- loglik(c(x[i], y[j]), sigma = Sigma)
## }
## persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

