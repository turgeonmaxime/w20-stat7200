## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=FALSE,
                      message = FALSE)
library(tidyverse)


## ----echo = FALSE-------------------------------------------------------------
C <- chol(S <- toeplitz(.9 ^ (0:31))) # Cov.matrix and its root
set.seed(17)
X <- matrix(rnorm(32000), 1000, 32)
Z <- X %*% C  ## ==>  cov(Z) ~=  C'C = S
pZ <- prcomp(Z, tol = 0.1)
screeplot(pZ, type = 'lines', main = '')


## -----------------------------------------------------------------------------
library(mvtnorm)
Sigma <- matrix(c(1, 0.5, 0.1,
                  0.5, 1, 0.5,
                  0.1, 0.5, 1),
                ncol = 3) 

set.seed(17)
X <- rmvnorm(100, sigma = Sigma)
pca <- prcomp(X)


## -----------------------------------------------------------------------------
summary(pca)


## -----------------------------------------------------------------------------
screeplot(pca, type = 'l')


## -----------------------------------------------------------------------------
pca <- prcomp(USArrests, scale = TRUE)


## -----------------------------------------------------------------------------
summary(pca)


## -----------------------------------------------------------------------------
screeplot(pca, type = 'l')


## -----------------------------------------------------------------------------
library(mvtnorm)
Sigma <- matrix(c(1, 0.5, 0.1,
                  0.5, 1, 0.5,
                  0.1, 0.5, 1),
                ncol = 3) 
mu <- c(1, 2, 2)


## -----------------------------------------------------------------------------
set.seed(17)
X <- rmvnorm(100, mean = mu,
             sigma = Sigma)
pca <- prcomp(X)

colMeans(X)
colMeans(pca$x)


## -----------------------------------------------------------------------------
# On the other hand
pca <- prcomp(X, center = FALSE)
colMeans(pca$x)


## ----mnist, cache = TRUE, message=FALSE---------------------------------------
library(tidyverse)
library(dslabs)

mnist <- read_mnist()

dim(mnist$train$images)
dim(mnist$test$images)

head(mnist$train$labels)


## -----------------------------------------------------------------------------
img <- matrix(mnist$train$images[1,], ncol = 28) 
# Switch columns
img <- img[,rev(1:28)]
image(img, col = gray.colors(12, rev = TRUE),
      axes = FALSE, asp = 1)


## ----mnistPCA, cache = TRUE---------------------------------------------------
decomp <- prcomp(mnist$train$images)


## -----------------------------------------------------------------------------
screeplot(decomp, type = 'lines',
          npcs = 20, main = "")


## -----------------------------------------------------------------------------
decomp$x[,1:2] %>%
  as.data.frame() %>%
  mutate(label = factor(mnist$train$labels)) %>%
  ggplot(aes(PC1, PC2, colour = label)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  coord_fixed()


## -----------------------------------------------------------------------------
# And on the test set
decomp %>%
  predict(newdata = mnist$test$images) %>%
  as.data.frame() %>%
  mutate(label = factor(mnist$test$labels)) %>%
  ggplot(aes(PC1, PC2, colour = label)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  coord_fixed()


## -----------------------------------------------------------------------------
par(mfrow = c(2, 2))
for (i in seq_len(4)) {
  img <- matrix(decomp$rotation[,i], ncol = 28)
  img <- img[,rev(1:28)]
  image(img, col = gray.colors(12, rev = TRUE),
        axes = FALSE, main = paste0("PC", i), asp = 1)
}


## -----------------------------------------------------------------------------
# Approximation with 90 PCs
approx_mnist <- decomp$rotation[, seq_len(90)] %*% 
  decomp$x[1, seq_len(90)]

# Original image
img1 <- matrix(mnist$train$images[1,], ncol = 28)
img1 <- img1[,rev(1:28)]
# Approximation
img2 <- matrix(approx_mnist, ncol = 28)
img2 <- img2[,rev(1:28)]


## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
image(img1, col = gray.colors(12, rev = TRUE),
      axes = FALSE, main = "Original", asp = 1)
image(img2, col = gray.colors(12, rev = TRUE),
      axes = FALSE, main = "Approx", asp = 1)


## -----------------------------------------------------------------------------
set.seed(1234)
# Random measurement error
sigma <- 5

# Exact relationship between 
# Celsius and Fahrenheit
temp_c <- seq(-40, 40, by = 1)
temp_f <- 1.8*temp_c + 32


## -----------------------------------------------------------------------------
# Add measurement error
temp_c_noise <- temp_c + rnorm(n = length(temp_c), 
                               sd = sigma)
temp_f_noise <- temp_f + rnorm(n = length(temp_f), 
                               sd = sigma)


## -----------------------------------------------------------------------------
# Linear model
(fit <- lm(temp_f_noise ~ temp_c_noise))
confint(fit)

# PCA
pca <- prcomp(cbind(temp_c_noise, temp_f_noise))
pca$rotation
pca$rotation[2,"PC1"]/pca$rotation[1,"PC1"]


## ---- message = FALSE---------------------------------------------------------
library(dslabs)
library(ggridges)

# Data on Breast Cancer
as.data.frame(brca$x) %>% 
  gather(variable, measurement) %>% 
  mutate(variable = reorder(variable, measurement, 
                            median)) %>% 
  ggplot(aes(x = measurement, y = variable)) + 
  geom_density_ridges() + theme_ridges() +
  coord_cartesian(xlim = c(0, 250))


## -----------------------------------------------------------------------------
# Remove some variables
rem_index <- which(colnames(brca$x) %in%
                     c("area_worst", "area_mean",
                       "perimeter_worst", 
                       "perimeter_mean"))
dataset <- brca$x[,-rem_index]
decomp <- prcomp(dataset)
summary(decomp)$importance[,1:3]


## -----------------------------------------------------------------------------
screeplot(decomp, type = 'l')


## -----------------------------------------------------------------------------
# Let's put a CI around the first eigenvalue
first_ev <- decomp$sdev[1]^2
n <- nrow(dataset)

# Recall that TV = 2166
c("LB" = first_ev/(1+qnorm(0.975)*sqrt(2/n)),
  "Est." = first_ev,
  "UP" = first_ev/(1-qnorm(0.975)*sqrt(2/n)))


## -----------------------------------------------------------------------------
B <- 1000; n <- 100; p <- 3

results <- purrr::map_df(seq_len(B), function(b) {
    X <- matrix(rnorm(p*n, sd = sqrt(c(1, 2, 3))), 
                ncol = p, byrow = TRUE)
    tmp <- eigen(cov(X), symmetric = TRUE, 
                 only.values = TRUE)
    tibble(ev1 = tmp$values[1],
           ev2 = tmp$values[2],
           ev3 = tmp$values[3])
})


## -----------------------------------------------------------------------------
results %>% 
  gather(ev, value) %>% 
  ggplot(aes(value, fill = ev)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  geom_vline(xintercept = c(1, 2, 3),
             linetype = 'dashed')


## -----------------------------------------------------------------------------
summarise_all(results, mean)


## -----------------------------------------------------------------------------
# Is there some correlation?
GGally::ggpairs(results)


## -----------------------------------------------------------------------------
p <- 2
results_vect <- purrr::map_df(seq_len(B), function(b) {
  X <- matrix(rnorm(p*n, sd = c(1, 2)), ncol = p,
              byrow = TRUE)
  tmp <- eigen(cov(X), symmetric = TRUE)
  
  tibble(
    xend = tmp$vectors[1,1],
    yend = tmp$vectors[2,1]
  )
})


## -----------------------------------------------------------------------------
ggplot(results_vect) +
  geom_segment(aes(xend = xend, yend = yend),
               x = 0, y = 0, colour = 'grey60') +
  geom_segment(x = 0, xend = 0,
               y = 0, yend = 1,
               colour = 'blue', size = 2) +
  expand_limits(y = 0, x = c(-1, 1)) +
  theme_minimal()


## ---- message = FALSE---------------------------------------------------------
# Or looking at angles
transmute(results_vect, 
          theta = atan2(yend, xend)) %>% 
  ggplot(aes(theta)) +
  geom_histogram() +
  theme_minimal() + 
  geom_vline(xintercept = pi/2)


## -----------------------------------------------------------------------------
library(mvtnorm)
# What about the t distribution
results_t <- purrr::map_df(seq_len(B), function(b) {
    X <- rmvt(n, sigma = diag(2*c(1, 2, 3)/4), 
              df = 4)
    tmp <- eigen(cov(X), symmetric = TRUE, 
                 only.values = TRUE)
    tibble(ev1 = tmp$values[1],
           ev2 = tmp$values[2],
           ev3 = tmp$values[3])
})


## -----------------------------------------------------------------------------
results_t %>% 
  gather(ev, value) %>% 
  ggplot(aes(value, fill = ev)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  geom_vline(xintercept = c(1, 2, 3),
             linetype = 'dashed')


## -----------------------------------------------------------------------------
summarise_all(results_t, mean)


## ---- message = FALSE, echo = FALSE-------------------------------------------
library(cowplot)
plot1 <- results %>% 
  gather(ev, value) %>% 
  ggplot(aes(value, fill = ev)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = 'none') +
  geom_vline(xintercept = c(1, 2, 3),
             linetype = 'dashed') +
  coord_cartesian(xlim = c(0, 10),
                  ylim = c(0, 3)) +
  ggtitle("Normal dist.") + ylab("")

plot2 <- results_t %>% 
  gather(ev, value) %>% 
  ggplot(aes(value, fill = ev)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = 'none') +
  geom_vline(xintercept = c(1, 2, 3),
             linetype = 'dashed') +
  coord_cartesian(xlim = c(0, 10),
                  ylim = c(0, 3)) +
  ggtitle("t dist. (df = 4)") + ylab("")

cowplot::plot_grid(plot1, plot2)


## -----------------------------------------------------------------------------
# Keep only mean measurements
rem_index2 <- grep("mean", colnames(brca$x),
                   invert = TRUE)
dataset <- brca$x[,-c(rem_index,
                      rem_index2)]
R <- cor(dataset)

pairs(dataset)


## -----------------------------------------------------------------------------
# Overall mean
r_bar <- mean(R[upper.tri(R, diag = FALSE)])

# Column specific means
r_cols <- (colSums(R) - 1)/(nrow(R) - 1)

# Extra quantities
p <- ncol(dataset)
n <- nrow(dataset)
gamma_hat <- (p - 1)^2*(1 - (1 - r_bar)^2)/
  (p - (p - 2)*(1 - r_bar)^2)


## -----------------------------------------------------------------------------
# Test statistic
Tstat <- sum((R[upper.tri(R, 
                          diag = FALSE)] - r_bar)^2) - 
  gamma_hat*sum((r_cols - r_bar)^2)
Tstat <- (n-1)*Tstat/(1-r_bar)^2

Tstat > qchisq(0.95, 0.5*(p+1)*(p-2))


## -----------------------------------------------------------------------------
decomp <- prcomp(dataset)

summary(decomp)$importance[,seq_len(3)]
screeplot(decomp, type = 'l')


## -----------------------------------------------------------------------------
permute_data <- function(data) {
  p <- ncol(data)
  data_perm <- data
  for (i in seq_len(p)) {
    ind_sc <- sample(nrow(data))
    data_perm[,i] <- data[ind_sc, i]
  }
  return(data_perm)
}


## -----------------------------------------------------------------------------
set.seed(123)
B <- 1000
alpha <- 0.05
results <- matrix(NA, ncol = B, 
                  nrow = ncol(dataset))

results[,1] <- decomp$sdev
results[,-1] <- replicate(B - 1, {
  data_perm <- permute_data(dataset)
  prcomp(data_perm)$sdev
})


## -----------------------------------------------------------------------------
cutoff <- apply(results, 1, function(row) {
  mean(row >= row[1])
})
which(cutoff < alpha)


## -----------------------------------------------------------------------------
# Continuing with our example on breast cancer
decomp <- prcomp(dataset)

# Extract PCs and loadings
PCs <- decomp$x[, 1:2]
loadings <- decomp$rotation[, 1:2]

# Extract data on tumour type
colour <- ifelse(brca$y == "B", "black", 'blue')


## -----------------------------------------------------------------------------
par(mfrow = c(1,2))
plot(PCs, pch = 19, col = colour)
plot(loadings, type = 'n')
text(loadings, 
     labels = colnames(dataset),
     col = 'red')
arrows(0, 0, 0.9 * loadings[, 1], 
       0.9 * loadings[, 2], 
       col = 'red', 
       length = 0.1)


## ---- echo = -1, message = FALSE, warning = FALSE-----------------------------
par(mfrow = c(1, 1))
# Or both on the same plot
plot(PCs, pch = 19, col = colour)
text(loadings, 
     labels = colnames(dataset),
     col = 'red')
arrows(0, 0, 0.9 * loadings[, 1], 
       0.9 * loadings[, 2], 
       col = 'red', 
       length = 0.1)


## -----------------------------------------------------------------------------
# The biplot function rescales for us
biplot(decomp)


## -----------------------------------------------------------------------------
# With scaled data
biplot(prcomp(dataset, scale = TRUE))


## ----message = FALSE----------------------------------------------------------
library(tidyverse)
url <- "https://maxturgeon.ca/w20-stat7200/prostate.csv"
prostate <- read_csv(url)

# Separate into training and testing sets
data_train <- filter(prostate, train == TRUE) %>%
  dplyr::select(-train)
data_test <- filter(prostate, train == FALSE) %>%
  dplyr::select(-train)


## ----message = FALSE----------------------------------------------------------
# First model: Linear regression
lr_model <- lm(lpsa ~ ., data = data_train)
lr_pred <- predict(lr_model, newdata = data_test)
(lr_mse <- mean((data_test$lpsa - lr_pred)^2))

# PCA
decomp <- data_train %>%
  subset(select = -lpsa) %>%
  as.matrix() %>%
  prcomp
summary(decomp)$importance[,1:3]

screeplot(decomp, type = 'lines')


## -----------------------------------------------------------------------------
# Second model: PCs for predictors
train_pc <- data_train
train_pc$PC1 <- decomp$x[,1]
pc_model <- lm(lpsa ~ PC1, data = train_pc)


## -----------------------------------------------------------------------------
test_pc <- as.data.frame(predict(decomp, data_test))
pc_pred <- predict(pc_model,
                   newdata = test_pc)
(pc_mse <- mean((data_test$lpsa - pc_pred)^2))


## -----------------------------------------------------------------------------
contribution <- decomp$rotation[,"PC1"]
round(contribution, 3)[1:6]
round(contribution, 3)[7:8]

(keep <- names(which(abs(contribution) > 0.01)))

fs_model <- lm(lpsa ~ ., data = data_train[,c(keep, "lpsa")])
fs_pred <- predict(fs_model, newdata = data_test)
(fs_mse <- mean((data_test$lpsa - fs_pred)^2))


## ---- message = FALSE---------------------------------------------------------
model_plot <- data.frame(
  "obs" = data_test$lpsa,
  "LR" = lr_pred,
  "PC" = pc_pred,
  "FS" = fs_pred
) %>%
  gather(Model, pred, -obs)


## ---- message = FALSE---------------------------------------------------------
ggplot(model_plot,
       aes(pred, obs, colour = Model)) +
  geom_point() +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0) +
  theme(legend.position = 'top') +
  xlab("Predicted") + ylab("Observed")

