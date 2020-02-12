## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=FALSE)
options(knitr.kable.NA = '-')
library(pander)

panderOptions('round', 2)
panderOptions('keep.trailing.zeros', TRUE)
# No caption
set.caption("", permanent = TRUE)


## ---- message = FALSE---------------------------------------------------------
# Let's revisit the plastic film data
library(heplots)
library(tidyverse)

Y <- Plastic %>% 
  select(tear, gloss, opacity) %>% 
  as.matrix

X <- model.matrix(~ rate, data = Plastic)
head(X)


## -----------------------------------------------------------------------------
(B_hat <- solve(crossprod(X)) %*% t(X) %*% Y)


## -----------------------------------------------------------------------------
# Compare with lm output
fit <- lm(cbind(tear, gloss, opacity) ~ rate, 
          data = Plastic)
coef(fit)


## -----------------------------------------------------------------------------
Y_hat <- fitted(fit)
residuals <- residuals(fit)

crossprod(Y_hat, residuals)
crossprod(X, residuals)

# Is this really zero?
isZero <- function(mat) {
  all.equal(mat, matrix(0, ncol = ncol(mat), 
                        nrow = nrow(mat)),
            check.attributes = FALSE)
}

isZero(crossprod(Y_hat, residuals))
isZero(crossprod(X, residuals))


## ---- message = FALSE---------------------------------------------------------
library(heplots)

head(NLSY)


## -----------------------------------------------------------------------------
# Fit model and look at coefficients
fit <- lm(cbind(math, read) ~ income + educ, 
          data = NLSY)

coef(fit)


## -----------------------------------------------------------------------------
range(NLSY$income)
range(NLSY$educ)


## -----------------------------------------------------------------------------
# Recall our model for Plastic
fit <- lm(cbind(tear, gloss, opacity) ~ rate, 
          data = Plastic)

new_x <- data.frame(rate = factor("High", 
                                  levels = c("Low", 
                                             "High")))
(prediction <- predict(fit, newdata = new_x))


## -----------------------------------------------------------------------------
X <- model.matrix(fit)
S <- crossprod(resid(fit))/(nrow(Plastic) - ncol(X))
new_x <- model.matrix(~rate, new_x)

quad_form <- drop(new_x %*% solve(crossprod(X)) %*% 
                    t(new_x))

# Estimation covariance
(est_cov <- S * quad_form) 

# Forecasting covariance
(fct_cov <- S *(1 + quad_form)) 


## -----------------------------------------------------------------------------
# Estimation CIs
cbind(drop(prediction) - 1.96*sqrt(diag(est_cov)),
      drop(prediction) + 1.96*sqrt(diag(est_cov)))

# Forecasting CIs
cbind(drop(prediction) - 1.96*sqrt(diag(fct_cov)),
      drop(prediction) + 1.96*sqrt(diag(fct_cov)))


## ---- warning = FALSE, digits = 2---------------------------------------------
# Going back to our NLSY example
full_model <- lm(cbind(math, read) ~ income + educ + 
                   antisoc + hyperact,
                 data = NLSY)

library(pander)
pander(anova(full_model, test = "Wilks"))

pander(anova(full_model, test = "Roy"))


## -----------------------------------------------------------------------------
# Visualize the error and hypothesis ellipses
heplot(full_model)


## ---- warning = FALSE---------------------------------------------------------
# Fit a model with only income and educ
rest_model <- lm(cbind(math, read) ~ income + educ,
                 data = NLSY)

pander(anova(full_model, rest_model,
             test = "Wilks"))

pander(anova(full_model, rest_model,
             test = "Roy"))


## -----------------------------------------------------------------------------
# Let's look at the eigenvalues
E <- crossprod(residuals(full_model))
H <- crossprod(residuals(rest_model)) - E

result <- eigen(H %*% solve(E),
                only.values = TRUE)
result$values[seq_len(2)]


## ----eval = -1----------------------------------------------------------------
AIC(full_model)
# Error in logLik.lm(full_model) : 
#   'logLik.lm' does not support multiple responses
class(full_model)


## -----------------------------------------------------------------------------
logLik.mlm <- function(object, ...) {
  resids <- residuals(object)
  Sigma_ML <- crossprod(resids)/nrow(resids)
  ans <- sum(mvtnorm::dmvnorm(resids, log = TRUE,
                              sigma = Sigma_ML))
  df <- prod(dim(coef(object))) + 
    choose(ncol(Sigma_ML) + 1, 2)
  attr(ans, "df") <- df
  class(ans) <- "logLik"
  return(ans)
}


## -----------------------------------------------------------------------------
logLik(full_model)

AIC(full_model)
AIC(rest_model)


## -----------------------------------------------------------------------------
# Model selection for Plastic data
lhs <- "cbind(tear, gloss, opacity) ~"
rhs_form <- c("1", "rate", "additive", 
              "rate+additive", "rate*additive")

purrr::map_df(rhs_form, function(rhs) {
  form <- formula(paste(lhs, rhs))
  fit <- lm(form, data = Plastic)
  return(data.frame(model = rhs, aic = AIC(fit),
                    stringsAsFactors = FALSE))
})


## ---- message = FALSE, echo = FALSE, eval = FALSE-----------------------------
## # Full subset selection for NLSY data
## library(leaps)
## # DOESNT WORK
## out <- regsubsets(cbind(math, read) ~ income + educ +
##                    antisoc + hyperact,
##                  data = NLSY,
##                  nbest = 1,
##                  nvmax = NULL,
##                  force.in = NULL, force.out = NULL,
##                  method = "exhaustive")


## ---- message = FALSE---------------------------------------------------------
library(openintro)
model <- lm(cbind(startPr, totalPr) ~ 
              nBids + cond + sellerRate + 
              wheels + stockPhoto, 
            data = marioKart)

X <- model.matrix(model)
P <- X %*% solve(crossprod(X)) %*% t(X)
lev_values <- diag(P)

hist(lev_values, 50)


## -----------------------------------------------------------------------------
n <- nrow(marioKart)
resids <- residuals(model)
S <- crossprod(resids)/(n - ncol(X))

S_inv <- solve(S)

const <- lev_values/((1 - lev_values)^2*ncol(X))
cook_values <- const * diag(resids %*% S_inv 
                            %*% t(resids))

hist(cook_values, 50)


## -----------------------------------------------------------------------------
# Cut-off value
(cutoff <- qchisq(0.5, ncol(S)*(n - ncol(X))))
which(cook_values > cutoff)

