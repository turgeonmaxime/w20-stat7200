library(heplots)
library(tidyverse)

# Multivariate Likelihood----
# We need to define this in order for R
# to automatically compute the AIC
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

# Fit full and empty models----
full_model <- lm(cbind(math, read) ~ income + educ + 
                   antisoc + hyperact,
                 data = NLSY)
null_model <- lm(cbind(math, read) ~ 1,
                 data = NLSY)

# Compute their AIC
AIC(full_model)
AIC(null_model)

# Generate all postible models----
list_covs <- c("income", "educ", "antisoc", "hyperact")
lhs <- "cbind(math, read)"

# Create formulas for each model
p <- length(list_covs)
list_models <- purrr::map(seq_len(p), function(k) {
  # Find all length k combinations of covariates
  list_fixed_k <- combn(list_covs, k, simplify = FALSE)
  # Go through combinations to create formulas
  purrr::map(list_fixed_k, function(covs) {
    rhs <- paste(covs, collapse = "+")
    paste(lhs, rhs, sep = "~")
  })
})

# Turn into vector
list_models <- unlist(list_models)

# Fit all possible models
fit_models <- purrr::map(list_models, function(form) {
  lm(as.formula(form),
     data = NLSY)
})

# Compute their AIC
AIC_models <- purrr::map_dbl(fit_models, AIC)

AIC_data <- data.frame(formula = list_models,
                       AIC = AIC_models,
                       stringsAsFactors = FALSE)

# Add null model
AIC_data <- bind_rows(
  AIC_data,
  data.frame(formula = paste(lhs, "1", sep = "~"),
             AIC = AIC(null_model),
             stringsAsFactors = FALSE)
)

# Find minimum AIC
dplyr::filter(AIC_data, AIC == min(AIC))
