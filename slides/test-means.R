## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=FALSE)


## ---- message=FALSE-----------------------------------------------------------
library(dslabs)
library(tidyverse)

dataset <- filter(gapminder, year == 2012, 
                  !is.na(infant_mortality))

dataset <- dataset[,c("infant_mortality",
                      "life_expectancy",
                      "fertility")]
dataset <- as.matrix(dataset)


## ---- message=FALSE-----------------------------------------------------------
dim(dataset)


## ---- message=FALSE-----------------------------------------------------------
# Assume we know Sigma
Sigma <- matrix(c(555, -170, 30, -170, 65, -10, 
                  30, -10, 2), ncol = 3)

mu_hat <- colMeans(dataset) 
mu_hat


## ---- message=FALSE-----------------------------------------------------------
# Test mu = mu_0
mu_0 <- c(25, 50, 3)
test_statistic <- nrow(dataset) * t(mu_hat - mu_0) %*% 
  solve(Sigma) %*% (mu_hat - mu_0)

c(drop(test_statistic), qchisq(0.95, df = 3))

drop(test_statistic) > qchisq(0.95, df = 3)


## ---- message=FALSE-----------------------------------------------------------
n <- nrow(dataset); p <- ncol(dataset)

# Test mu = mu_0
mu_0 <- c(25, 50, 3)
test_statistic <- n * t(mu_hat - mu_0) %*% 
  solve(cov(dataset)) %*% (mu_hat - mu_0)

critical_val <- (n - 1)*p*qf(0.95, df1 = p,
                             df2 = n - p)/(n-p)


## ---- message=FALSE-----------------------------------------------------------
c(drop(test_statistic), critical_val)

drop(test_statistic) > critical_val


## -----------------------------------------------------------------------------
n <- nrow(dataset); p <- ncol(dataset)

critical_val <- (n - 1)*p*qf(0.95, df1 = p,
                             df2 = n - p)/(n-p)
sample_cov <- diag(cov(dataset))


## -----------------------------------------------------------------------------
cbind(mu_hat - sqrt(critical_val*
                      sample_cov/n),
      mu_hat + sqrt(critical_val*
                      sample_cov/n))


## -----------------------------------------------------------------------------
U <- matrix(c(1, 0, 0,
              0, 1, 0), 
            ncol = 2)
R <- n*solve(t(U) %*% cov(dataset) %*% U)
transf <- chol(R)


## -----------------------------------------------------------------------------
# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), 
                                     sin(theta_vect))
# Then turn into ellipse
ellipse <- circle %*% t(solve(transf)) + 
  matrix(mu_hat[1:2], ncol = 2, 
         nrow = nrow(circle), 
         byrow = TRUE)


## -----------------------------------------------------------------------------
# Eigendecomposition
# To visualize the principal axes
decomp <- eigen(t(U) %*% cov(dataset) %*% U)
first <- sqrt(decomp$values[1]) *
  decomp$vectors[,1] * sqrt(critical_val)
second <- sqrt(decomp$values[2]) * 
  decomp$vectors[,2] * sqrt(critical_val)


## ---- echo = FALSE------------------------------------------------------------
plot(ellipse, type = 'l', asp = 1,
     xlab = colnames(dataset)[1],
     ylab = colnames(dataset)[2])
lines(x = c(0 + mu_hat[1], first[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], first[2]/sqrt(n) + mu_hat[2]))
lines(x = c(0 + mu_hat[1], -first[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], -first[2]/sqrt(n) + mu_hat[2]))
lines(x = c(0 + mu_hat[1], second[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], second[2]/sqrt(n) + mu_hat[2]))
lines(x = c(0 + mu_hat[1], -second[1]/sqrt(n) + mu_hat[1]),
      y = c(0 + mu_hat[2], -second[2]/sqrt(n) + mu_hat[2]))
points(x = mu_hat[1],
       y = mu_hat[2])

# Add 1d projections
axis_proj <- cbind(mu_hat - sqrt(critical_val*
                                   sample_cov/n),
                   mu_hat + sqrt(critical_val*
                                   sample_cov/n))
abline(v = axis_proj[1,], lty = 2)
abline(h = axis_proj[2,], lty = 2)


## ---- message=FALSE-----------------------------------------------------------
# Let's focus on only two variables
dataset <- dataset[,c("infant_mortality",
                      "life_expectancy")]

n <- nrow(dataset); p <- ncol(dataset)


## -----------------------------------------------------------------------------
alpha <- 0.05
mu_hat <- colMeans(dataset) 
sample_cov <- diag(cov(dataset))

# Simultaneous CIs
critical_val <- (n - 1)*p*qf(1-0.5*alpha, df1 = p,
                             df2 = n - p)/(n-p)

simul_ci <- cbind(mu_hat - sqrt(critical_val*
                                  sample_cov/n),
                  mu_hat + sqrt(critical_val*
                                  sample_cov/n))


## -----------------------------------------------------------------------------
# Univariate without correction
univ_ci <- cbind(mu_hat - qt(1-0.5*alpha, n - 1) *
                   sqrt(sample_cov/n),
                 mu_hat + qt(1-0.5*alpha, n - 1) *
                   sqrt(sample_cov/n))


## -----------------------------------------------------------------------------
# Bonferroni adjustment
bonf_ci <- cbind(mu_hat - qt(1-0.5*alpha/p, n - 1) *
                   sqrt(sample_cov/n),
                 mu_hat + qt(1-0.5*alpha/p, n - 1) *
                   sqrt(sample_cov/n))


## -----------------------------------------------------------------------------
simul_ci
univ_ci
bonf_ci


## ---- echo = FALSE------------------------------------------------------------
transf_mat <- solve(chol(solve(cov(dataset)/n)))
# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
ellipse1 <- circle %*% t(transf_mat)
ellipse2 <- t(apply(ellipse1, 1, function(row) row + mu_hat))
colnames(ellipse2) <- c("X", "Y")

data_ellipse <- data.frame(ellipse2)

bind_rows(
  data.frame(t(simul_ci)) %>% mutate(type = 'T2-intervals'),
  data.frame(t(univ_ci)) %>% mutate(type = 'Unadjusted'),
  data.frame(t(bonf_ci)) %>% mutate(type = 'Bonferroni')
  ) %>% 
  ggplot() +
  geom_polygon(data = data_ellipse,
               aes(X, Y),
               fill = 'grey60')+
  geom_vline(aes(xintercept = infant_mortality,
                 linetype = type)) +
  geom_hline(aes(yintercept = life_expectancy,
                 linetype = type)) +
  theme_minimal() +
  geom_point(x = mu_hat[1],
             y = mu_hat[2],
             size = 2) +
  xlab("Infant mortality") +
  ylab("Life Expectancy") +
  theme(legend.position = "top",
        legend.title=element_blank()) +
  scale_linetype_discrete(breaks = c('T2-intervals',
                                     'Bonferroni',
                                     'Unadjusted')) +
  coord_fixed()


## ----message = FALSE----------------------------------------------------------
dataset1 <- filter(gapminder, year == 2012, 
                   continent == "Africa",
                   !is.na(infant_mortality))

dataset1 <- dataset1[,c("life_expectancy",
                        "infant_mortality")]
dataset1 <- as.matrix(dataset1)
dim(dataset1)

dataset2 <- filter(gapminder, year == 2012, 
                   continent == "Asia",
                   !is.na(infant_mortality))

dataset2 <- dataset2[,c("life_expectancy",
                        "infant_mortality")]
dataset2 <- as.matrix(dataset2)
dim(dataset2)

n1 <- nrow(dataset1); n2 <- nrow(dataset2)
p <- ncol(dataset1)


## -----------------------------------------------------------------------------
(mu_hat1 <- colMeans(dataset1))
(mu_hat2 <- colMeans(dataset2))

(S1 <- cov(dataset1))
(S2 <- cov(dataset2))

# Even though it doesn't look reasonable
# We will assume equal covariance for now


## -----------------------------------------------------------------------------
mu_hat_diff <- mu_hat1 - mu_hat2

S_pool <- ((n1 - 1)*S1 + (n2 - 1)*S2)/(n1+n2-2)

test_statistic <- t(mu_hat_diff) %*% 
  solve((n1^-1 + n2^-1)*S_pool) %*% mu_hat_diff

const <- (n1 + n2 - 2)*p/(n1 + n2 - p - 2)
critical_val <- const * qf(0.95, df1 = p,
                           df2 = n1 + n2 - p - 2)


## -----------------------------------------------------------------------------
c(drop(test_statistic), critical_val)

drop(test_statistic) > critical_val


## ----echo = FALSE-------------------------------------------------------------
R <- solve((n1^-1 + n2^-1)*S_pool)
transf <- chol(R)

# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse1 <- circle %*% t(solve(transf)) + 
  matrix(mu_hat_diff, ncol = p, 
         nrow = nrow(circle), 
         byrow = TRUE)

plot(ellipse1, type = 'l', asp = 1,
     xlab = colnames(dataset1)[1],
     ylab = colnames(dataset1)[2],
     main = "Comparing Africa vs. Asia")
points(x = mu_hat_diff[1],
       y = mu_hat_diff[2])


## -----------------------------------------------------------------------------
test_statistic <- t(mu_hat_diff) %*% 
  solve(n1^-1*S1 + n2^-1*S2) %*% mu_hat_diff

critical_val <- qchisq(0.95, df = p)

c(drop(test_statistic), critical_val)
drop(test_statistic) > critical_val


## ----echo = FALSE-------------------------------------------------------------
R <- solve(n1^-1*S1 + n2^-1*S2)
transf <- chol(R)

# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse2 <- circle %*% t(solve(transf)) + 
  matrix(mu_hat_diff, ncol = p, 
         nrow = nrow(circle), 
         byrow = TRUE)


## -----------------------------------------------------------------------------
W1 <- S1 %*% solve(n1^-1*S1 + n2^-1*S2)/n1
W2 <- S2 %*% solve(n1^-1*S1 + n2^-1*S2)/n2

trace_square <- sum(diag(W1%*%W1))/n1 + 
  sum(diag(W2%*%W2))/n2
square_trace <- sum(diag(W1))^2/n1 + 
  sum(diag(W2))^2/n2

(nu <- (p + p^2)/(trace_square + square_trace))


## -----------------------------------------------------------------------------
const <- nu*p/(nu - p - 1)
critical_val <- const * qf(0.95, df1 = p,
                           df2 = nu - p - 1)

c(drop(test_statistic), critical_val)
drop(test_statistic) > critical_val


## ----echo = FALSE-------------------------------------------------------------
# First create a circle of radius c
theta_vect <- seq(0, 2*pi, length.out = 100)
circle <- sqrt(critical_val) * cbind(cos(theta_vect), sin(theta_vect))
# Then turn into ellipse
ellipse3 <- circle %*% t(solve(transf)) + 
  matrix(mu_hat_diff, ncol = p, 
         nrow = nrow(circle), 
         byrow = TRUE)

xlim <- range(c(ellipse1[,1], 
                ellipse2[,1],
                ellipse3[,1]))
ylim <- range(c(ellipse1[,2], 
                ellipse2[,2],
                ellipse3[,2]))

plot(ellipse2, type = 'l', asp = 1,
     xlab = colnames(dataset1)[1],
     ylab = colnames(dataset1)[2],
     xlim = xlim, ylim = ylim,
     main = "Comparing Africa vs. Asia")
points(x = mu_hat_diff[1],
       y = mu_hat_diff[2])
lines(ellipse1, lty = 2)
lines(ellipse3, lty = 3)
legend('topright', legend = c("Unequal", "Equal", "Nel-VDM"), lty = 1:3)


## -----------------------------------------------------------------------------
library(mvtnorm)
set.seed(7200)

n <- 50; p <- 10
B <- 1000

# Simulate under the null
mu <- mu_0 <- rep(0, p)
# Cov: diag = 1; off-diag = 0.5
Sigma <- matrix(0.5, ncol = p, nrow = p)
diag(Sigma) <- 1


## -----------------------------------------------------------------------------
critical_val <- (n - 1)*p*qf(0.95, df1 = p,
                             df2 = n - p)/(n-p)

null_dist <- replicate(B, {
  Y_norm <- rmvnorm(n, mean = mu, sigma = Sigma)
  mu_hat <- colMeans(Y_norm)
  # Test mu = mu_0
  test_statistic <- n * t(mu_hat - mu_0) %*% 
    solve(cov(Y_norm)) %*% (mu_hat - mu_0)
})


## -----------------------------------------------------------------------------
# Type I error
mean(null_dist > critical_val)


## ----echo = FALSE-------------------------------------------------------------
# Actual density
x_vect <- 35*ppoints(B)
d_vect <- df(x_vect*(n-p)/((n - 1)*p), df1 = p,
             df2 = n - p)*(n-p)/((n - 1)*p)

hist(null_dist, 20, freq = FALSE,
     xlab = "Simulated data",
     main = "Black is smoothed density; Blue is theoretical density")
lines(density(null_dist), col = 'black')
lines(x_vect, d_vect, col = 'blue')


## -----------------------------------------------------------------------------
# Now the t distribution
nu <- 3

null_dist_t <- replicate(B, {
  Y_t <- rmvt(n, sigma = Sigma, df = nu, delta = mu)
  mu_hat <- colMeans(Y_t)
  # Test mu = mu_0
  test_statistic <- n * t(mu_hat - mu_0) %*% 
    solve(cov(Y_t)) %*% (mu_hat - mu_0)
})


## -----------------------------------------------------------------------------
# Type I error
mean(null_dist_t > critical_val)


## ----echo = FALSE-------------------------------------------------------------
# Actual density
hist(null_dist_t, 20, freq = FALSE,
     xlab = "Simulated data",
     main = "Black is smoothed density; Blue is theoretical density")
lines(density(null_dist_t), col = 'black')
lines(x_vect, d_vect, col = 'blue')


## -----------------------------------------------------------------------------
# Now a contaminated normal
sigma <- 3; epsilon <- 0.25
null_dist_cont <- replicate(B, {
  Z <- rmvnorm(n, sigma = diag(p))
  Y <- sample(c(sigma, 1), size = n, replace = TRUE,
              prob = c(epsilon, 1 - epsilon))*Z
  mu_hat <- colMeans(Y)
  # Test mu = mu_0
  test_statistic <- n * t(mu_hat - mu_0) %*% 
    solve(cov(Y)) %*% (mu_hat - mu_0)
})


## -----------------------------------------------------------------------------
# Type I error
mean(null_dist_cont > critical_val)


## ----echo = FALSE-------------------------------------------------------------
# Actual density
hist(null_dist_cont, 20, freq = FALSE,
     xlab = "Simulated data",
     main = "Black is smoothed density; Blue is theoretical density")
lines(density(null_dist_cont), col = 'black')
lines(x_vect, d_vect, col = 'blue')


## ----echo = FALSE, cache = TRUE-----------------------------------------------
library(tidyverse)
B <- 10000

# Simulate for several values of n
results <- purrr::map_df(seq(50, 300, by = 50),
                         function(n) {
  # First normal
  null_dist <- replicate(B, {
    Y_norm <- rmvnorm(n, mean = mu, sigma = Sigma)
    mu_hat <- colMeans(Y_norm)
    # Test mu = mu_0
    test_statistic <- n * t(mu_hat - mu_0) %*% 
      solve(cov(Y_norm)) %*% (mu_hat - mu_0)
  })
  
  # Second t distribution
  null_dist_t <- replicate(B, {
    Y_t <- rmvt(n, sigma = Sigma, df = nu, delta = mu)
    mu_hat <- colMeans(Y_t)
    # Test mu = mu_0
    test_statistic <- n * t(mu_hat - mu_0) %*% 
      solve(cov(Y_t)) %*% (mu_hat - mu_0)
  })
  
  # Third a contaminated normal
  null_dist_cont <- replicate(B, {
    Z <- rmvnorm(n, sigma = diag(p))
    Y <- sample(c(sigma, 1), size = n, replace = TRUE,
                prob = c(epsilon, 1 - epsilon))*Z
    mu_hat <- colMeans(Y)
    # Test mu = mu_0
    test_statistic <- n * t(mu_hat - mu_0) %*% 
      solve(cov(Y)) %*% (mu_hat - mu_0)
  })
  
  # Return in a data.frame
  bind_rows(
    mutate(tibble(statistics = null_dist), model = "Normal"),
    mutate(tibble(statistics = null_dist_t), model = "t-dist"),
    mutate(tibble(statistics = null_dist_cont), model = "Contaminated")
  ) %>% 
    mutate(n = n)
})


## ----echo = FALSE-------------------------------------------------------------
# Compute Type I error rate and plot as function of n
results %>% 
  group_by(model, n) %>% 
  # summarise(TypeI_fdist = mean(statistics > (n - 1)*p*qf(0.95, df1 = p,
  #                                                        df2 = n - p)/(n-p)),
  #           TypeI_chisq = mean(statistics > qchisq(0.95, df = p))) %>% 
  # gather(Test, TypeI, TypeI_fdist, TypeI_chisq) %>% 
  # mutate(Test = str_replace(Test, "TypeI_", "")) %>% 
  summarise(TypeI = mean(statistics > (n - 1)*p*qf(0.95, df1 = p,
                                                   df2 = n - p)/(n-p))) %>%
  ggplot(aes(n, TypeI, colour = model)) +
  geom_line() + geom_point() +
  # facet_grid(. ~ Test) %>% 
  theme_minimal() +
  theme(legend.position = 'top') +
  geom_hline(yintercept = 0.05, linetype = 'dashed') +
  expand_limits(y = 0) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))


## ---- cache = FALSE-----------------------------------------------------------
# Recall our dataset
dataset <- filter(gapminder, year == 2012, 
                  !is.na(infant_mortality))

dataset <- dataset[,c("infant_mortality",
                      "life_expectancy")]
dataset <- as.matrix(dataset)


## -----------------------------------------------------------------------------
# The sample estimators
colMeans(dataset)
cov(dataset)


## ---- message = FALSE---------------------------------------------------------
# The MCD estimators
library(rrcov)

mcd <- CovMcd(dataset)
getCenter(mcd)
getCov(mcd)


## -----------------------------------------------------------------------------
n <- nrow(dataset); p <- ncol(dataset)

# Classical T2
mu_0 <- c(25, 70)
test_statistic <- n * t(mu_hat - mu_0) %*% 
  solve(cov(dataset)) %*% (mu_hat - mu_0)

critical_val <- (n - 1)*p*qf(0.95, df1 = p,
                             df2 = n - p)/(n-p)


## -----------------------------------------------------------------------------
c(drop(test_statistic), critical_val)
drop(test_statistic) > critical_val


## -----------------------------------------------------------------------------
# Robust T2
t2_robust <- T2.test(dataset, mu = mu_0, method = "mcd")
t2_robust
t2_robust$p.value

