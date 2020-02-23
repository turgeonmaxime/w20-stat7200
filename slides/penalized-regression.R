## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=FALSE)


## ---- message = FALSE---------------------------------------------------------
library(tidyverse)
url <- "https://maxturgeon.ca/w20-stat7200/prostate.csv"
prostate <- read_csv(url)

# Separate into training and testing sets
data_train <- filter(prostate, train == TRUE) %>% 
  dplyr::select(-train)
data_test <- filter(prostate, train == FALSE) %>% 
  dplyr::select(-train)


## -----------------------------------------------------------------------------
# OLS
model1 <- lm(lpsa ~ ., 
             data = data_train)
pred1 <- predict(model1, data_test)

mean((data_test$lpsa - pred1)^2)


## -----------------------------------------------------------------------------
# Ridge regression
X_train <- model.matrix(lpsa ~ ., 
             data = data_train)
Y_train <- data_train$lpsa

B_ridge <- solve(crossprod(X_train) + diag(0.7, 9), 
                 t(X_train)) %*% Y_train


## -----------------------------------------------------------------------------
X_test <- model.matrix(lpsa ~ ., 
             data = data_test)

pred2 <- X_test %*% B_ridge

mean((data_test$lpsa - pred2)^2)


## -----------------------------------------------------------------------------
# Compare both estimates
head(cbind(coef(model1), B_ridge))


## -----------------------------------------------------------------------------
mse_df <- purrr::map_df(seq(0, 5, by = 0.1),
                        function(lambda) {
  B_ridge <- solve(crossprod(X_train) + diag(lambda, 9), 
                 t(X_train)) %*% Y_train
  pred2 <- X_test %*% B_ridge

  mse <- mean((data_test$lpsa - pred2)^2)
  return(data.frame(MSE = mse,
                    lambda = lambda))
  })


## -----------------------------------------------------------------------------
ols_mse <- mean((data_test$lpsa - pred1)^2)

ggplot(mse_df, aes(lambda, MSE)) + 
  geom_line() + theme_minimal() + 
  geom_hline(yintercept = ols_mse)


## ---- message = FALSE---------------------------------------------------------
library(glmnet)

# Fit for multiple values of lambda
X_train <- model.matrix(lpsa ~ . - 1,
                        data = data_train)
ridge_fit <- glmnet(X_train, data_train$lpsa,
                    alpha = 0,
                    lambda = seq(0, 5, by = 0.1))


## ---- message = FALSE---------------------------------------------------------
# Plot the value of the coefficients
# as a function of lambda
plot(ridge_fit, xvar = "lambda")
abline(h = 0, lty = 2)


## -----------------------------------------------------------------------------
# Fit lasso regression along the same lambda sequence
lasso_fit <- glmnet(X_train, data_train$lpsa, 
                    alpha = 1, # For lasso regression
                    lambda = seq(0, 5, by = 0.1))


## -----------------------------------------------------------------------------
X_test <- model.matrix(lpsa ~ . - 1,
                       data = data_test)
lasso_pred <- predict(lasso_fit, newx = X_test)
lasso_mse <- apply(lasso_pred, 2, function(col) {
  mean((data_test$lpsa - col)^2)
})


## -----------------------------------------------------------------------------
lasso_mse_df <- data.frame(MSE = lasso_mse,
                           lambda = seq(0, 5, by = 0.1))
ggplot(mse_df, aes(lambda, MSE)) + 
  geom_line() + theme_minimal() + 
  geom_hline(yintercept = ols_mse) +
  geom_line(data = lasso_mse_df, colour = 'red')


## -----------------------------------------------------------------------------
# Plot the value of the coefficients
# as a function of lambda
plot(lasso_fit, xvar = "lambda")
abline(h = 0, lty = 2)


## -----------------------------------------------------------------------------
# Where is the min MSE?
filter(lasso_mse_df, MSE == min(MSE))
# What are the estimates?
coef(lasso_fit, s = 4.9)


## -----------------------------------------------------------------------------
# Take all the data
dataset <- dplyr::select(prostate, -train)
dim(dataset)


## ---- message = FALSE---------------------------------------------------------
set.seed(7200)
library(caret)
# 5-fold CV
trainIndex <- createFolds(dataset$lpsa, k = 5)
str(trainIndex)


## -----------------------------------------------------------------------------
# Define function to compute MSE
compute_mse <- function(prediction, actual) {
  # Recall: the prediction comes in an array
  apply(prediction, 2, function(col) {
    mean((actual - col)^2)
    })
}


## -----------------------------------------------------------------------------
MSEs <- sapply(trainIndex, function(indices){
  X_train <- model.matrix(lpsa ~ . - 1,
                          data = dataset[-indices,])
  Y_train <- dataset$lpsa[-indices]
  X_test <- model.matrix(lpsa ~ . - 1,
                          data = dataset[indices,])
  lasso_fit <- glmnet(X_train, Y_train, alpha = 1, 
                      lambda = seq(0, 5, by = 0.1))
  lasso_pred <- predict(lasso_fit, newx = X_test)
  compute_mse(lasso_pred, dataset$lpsa[indices])
  })


## -----------------------------------------------------------------------------
# Each column is for a different fold
dim(MSEs)
CV_MSE <- colMeans(MSEs)

seq(0, 5, by = 0.1)[which.min(CV_MSE)]


## -----------------------------------------------------------------------------
# What are the estimates?
coef(lasso_fit, s = 0.4)


## -----------------------------------------------------------------------------
# Conveniently, glmnet has a function for CV
# It also chooses the lambda sequence for you
X <- model.matrix(lpsa ~ . -1, data = dataset)
lasso_cv_fit <- cv.glmnet(X, dataset$lpsa, alpha = 1,
                          nfolds = 5)


## -----------------------------------------------------------------------------
c("lambda.min" = lasso_cv_fit$lambda.min,
  "lambda.1se" = lasso_cv_fit$lambda.1se)

# What are the estimates?
coef(lasso_cv_fit, s = 'lambda.min')
# 1 SE rule
coef(lasso_cv_fit, s = 'lambda.1se')

