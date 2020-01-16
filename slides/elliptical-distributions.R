## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=TRUE, message = FALSE)


## -----------------------------------------------------------------------------
library(mvtnorm)

mu <- c(2, 1)
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
data <- expand.grid(seq(0, 4, length.out = 32),
                    seq(0, 2, length.out = 32))
data["dvalues"] <- dmvnorm(data, mean = mu, 
                           sigma = Sigma)


## -----------------------------------------------------------------------------
library(tidyverse)
ggplot(data, aes(Var1, Var2)) + 
  geom_contour(aes(z = dvalues))  +
  coord_fixed(ratio = 1)


## -----------------------------------------------------------------------------
k <- 0.12
const <- -2*log(k*2*pi*sqrt(det(Sigma)))

# Generate a circle
# First create a circle of radius const
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- const * cbind(cos(theta_vect), 
                        sin(theta_vect))


## -----------------------------------------------------------------------------
# Compute inverse Cholesky
transf_mat <- solve(chol(solve(Sigma)))
# Then turn circle into ellipse
ellipse <- circle %*% t(transf_mat)
# Then translate
ellipse <- t(apply(ellipse, 1, function(row) row + mu))


## -----------------------------------------------------------------------------
# Add ellipse to previous plot
ggplot(data, aes(Var1, Var2)) + 
  geom_contour(aes(z = dvalues)) +
  geom_polygon(data = data.frame(ellipse),
            aes(X1, X2), colour = 'red', fill = NA) +
  coord_fixed(ratio = 1)


## -----------------------------------------------------------------------------
set.seed(7200)

n <- 1000
p <- 2
Z <- rmvnorm(n, sigma = diag(p))


## -----------------------------------------------------------------------------
sigma <- 2
epsilon <- 0.25
w <- sample(c(sigma, 1), size = n, replace = TRUE,
            prob = c(epsilon, 1 - epsilon))

Y <- w*Z


## -----------------------------------------------------------------------------
# Plot the results
par(mfrow = c(1, 2))
plot(Z, main = 'Standard Normal',
     xlim = c(-4.5, 6.5), ylim = c(-7, 5))
plot(Y, main = 'Contaminated Normal',
     xlim = c(-4.5, 6.5), ylim = c(-7, 5))


## -----------------------------------------------------------------------------
# Colour points of Y according to 
# which distribution they come from
plot(Y, main = 'Contaminated Normal', col = w, 
     xlim = c(-4.5, 6.5), ylim = c(-7, 5))
legend("bottomleft", legend = c("Sigma = 1", "Sigma = 2"),
       col = c(1, 2), pch = 1)


## -----------------------------------------------------------------------------
# Let's look at the distribution of the sample means
B <- 1000; n <- 100
data <- purrr::map_df(seq_len(B), function(b) {
  Z <- rmvnorm(n, sigma = diag(p))
  Y <- sample(c(sigma, 1), size = n, replace = TRUE,
            prob = c(epsilon, 1 - epsilon)) * Z
  
  out <- data.frame(rbind(colMeans(Z), colMeans(Y)))
  out$Dist <- c("Standard", "Contaminated")
  return(out)
  })


## -----------------------------------------------------------------------------
ggplot(data, aes(X1, X2)) + 
  geom_point(aes(colour = Dist)) +
  theme(legend.position = 'top')


## -----------------------------------------------------------------------------
library(mvtnorm)
n <- 1000
mu <- c(1, 2)
Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
# Recall the multivariate case
Y_norm <- rmvnorm(n, mean = mu, sigma = Sigma)

colMeans(Y_norm)
cov(Y_norm)


## -----------------------------------------------------------------------------
# Now the t distribution
nu <- 4
Y_t <- rmvt(n, sigma = Sigma, df = nu, delta = mu)

colMeans(Y_t)
cov(Y_t)


## -----------------------------------------------------------------------------
data_plot <- rbind(
  mutate(data.frame(Y_norm), dist = "normal"),
  mutate(data.frame(Y_t), dist = "t")
)

ggplot(data_plot, aes(X1, X2)) +
  geom_point(alpha = 0.25) +
  geom_density_2d() +
  facet_grid(~ dist, labeller = label_both)

