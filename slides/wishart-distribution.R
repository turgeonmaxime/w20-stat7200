## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## -----------------------------------------------------------------------------
B <- 1000
n <- 10; p <- 4

traces <- replicate(B, {
  Z <- matrix(rnorm(n*p), ncol = p)
  W <- crossprod(Z)
  sum(diag(W))
})


## -----------------------------------------------------------------------------
hist(traces, 50, freq = FALSE)
lines(density(rchisq(B, df = n*p)))


## -----------------------------------------------------------------------------
B <- 1000
n <- 10
p <- 5

bartlett <- replicate(B, {
  X <- matrix(rnorm(n*p), ncol = p)
  L <- chol(crossprod(X))
})

dim(bartlett)


## ---- message = FALSE---------------------------------------------------------
library(tidyverse)

# Extract and plot diagonal^2
diagonal <- purrr::map_df(seq_len(B), function(i) {
  tmp <- diag(bartlett[,,i])^2
  data.frame(matrix(tmp, nrow = 1))
})


## -----------------------------------------------------------------------------
# Put into long format
diag_plot <- gather(diagonal, Entry, Value)

# Add chi-square means
diag_means <- data.frame(
  Entry = paste0("X", seq_len(p)),
  mean = n - seq_len(p) + 1
)


## -----------------------------------------------------------------------------
ggplot(diag_plot, aes(Value, fill = Entry)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  geom_vline(data = diag_means,
             aes(xintercept = mean, 
                 colour = Entry),
             linetype = 'dashed')


## ---- message = FALSE---------------------------------------------------------
# Extract and plot off-diagonal
off_diagonal <- purrr::map_df(seq_len(B), function(i) {
  tmp <- bartlett[,,i][upper.tri(bartlett[,,i])]
  data.frame(matrix(tmp, nrow = 1))
})
dim(off_diagonal)


## -----------------------------------------------------------------------------
# Put into long format
offdiag_plot <- gather(off_diagonal, Entry, Value)


## -----------------------------------------------------------------------------
ggplot(offdiag_plot, aes(Value, group = Entry)) +
  geom_density(fill = NA) +
  theme_minimal()


## ---- eval=TRUE, echo = FALSE-------------------------------------------------
# Ramus data, Timm (2002)
main_page <- "https://maxturgeon.ca/w20-stat7200/"
ramus <- read.csv(paste0(main_page, "Ramus.csv"))


## -----------------------------------------------------------------------------
var_names <- c("Age8", "Age8.5",
               "Age9", "Age9.5")

dataset <- ramus[,var_names]
dim(dataset)


## -----------------------------------------------------------------------------
# Sample covariance
Sn <- cov(dataset)

# Generalized variance
det(Sn)


## -----------------------------------------------------------------------------
# Simulate quantiles
set.seed(7200)
n <- nrow(dataset)
p <- ncol(dataset)
B <- 1000

simulated_vals <- replicate(B, {
  prod(rchisq(p, df = n - seq_len(p)))/((n-1)^p)
})


## -----------------------------------------------------------------------------
bounds <- quantile(simulated_vals,
                   probs = c(0.025, 0.975))

bounds


## -----------------------------------------------------------------------------
# 95% Confidence interval (reverse bounds)
det(Sn)/rev(bounds)


## -----------------------------------------------------------------------------
# Recall our covariance matrix for the Ramus dataset
round(Sn, 2)

# Visually we get
lattice::levelplot(Sn, xlab = "", ylab = "")


## -----------------------------------------------------------------------------
# Perhaps easier to interpret as correlations
# But be careful with the scale!
lattice::levelplot(cov2cor(Sn), 
                   xlab = "", ylab = "")



## -----------------------------------------------------------------------------
B <- 1000
n <- nrow(dataset)

boot_covs <- lapply(seq_len(B), function(b) {
  data_boot <- dataset[sample(n, n, replace = TRUE),]
  return(cov(data_boot))
})


## -----------------------------------------------------------------------------
# Extract the diagonal entries
diagonal <- purrr::map_df(boot_covs, function(Sn) {
  tmp <- diag(Sn)
  data.frame(matrix(tmp, nrow = 1))
  })


## -----------------------------------------------------------------------------
# Put into long format
diag_plot <- gather(diagonal, Entry, Value)

ggplot(diag_plot, aes(Value, fill = Entry)) +
  geom_density(alpha = 0.2) +
  theme_minimal()


## -----------------------------------------------------------------------------
# Multivariate normal theory predicts
# the diagonal entry should be scaled chi-square
ggplot(diag_plot, aes(sample = Value)) +
  geom_qq(distribution = qchisq, 
          dparams = list(df = n - 1)) +
  theme_minimal() + facet_wrap(~ Entry) + 
  geom_qq_line(distribution = qchisq, 
               dparams = list(df = n - 1))


## -----------------------------------------------------------------------------
# Finally, let's look at pairwise scatterplots 
# for off-diagonal entries
off_diag <- purrr::map_df(boot_covs, function(Sn) {
  tmp <- Sn[upper.tri(Sn)]
  data.frame(matrix(tmp, nrow = 1))
  })


## ---- message = FALSE---------------------------------------------------------
# Add column names
names(off_diag) <- c(paste0("8:",c("8.5","9","9.5")),
                     paste0("8.5:",c("9","9.5")), 
                     "9:9.5")

GGally::ggpairs(off_diag)

