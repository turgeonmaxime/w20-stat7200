## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## -----------------------------------------------------------------------------
library(copula)

# Gaussian copula where correlation is 0.5
gaus_copula <- normalCopula(0.5, dim = 2)
sample_copula1 <- rCopula(1000, gaus_copula)

plot(sample_copula1)

# Compare with independent copula, 
# i.e. two independent uniform variables.
gaus_copula <- normalCopula(0, dim = 2)
sample_copula2 <- rCopula(1000, gaus_copula)
plot(sample_copula2)


## ----echo = FALSE-------------------------------------------------------------
par(mfrow = c(1, 2))
plot(sample_copula1, main = "Corr. 0.5")
plot(sample_copula2, main = "Independent")


## -----------------------------------------------------------------------------
# Clayton copula with theta = 0.5
clay_copula <- claytonCopula(param = 0.5)
sample_copula1 <- rCopula(1000, clay_copula)

plot(sample_copula1)


## ----echo = FALSE, warning = FALSE, message = FALSE---------------------------
clay_copula <- claytonCopula(param = 0)
sample_copula1 <- rCopula(1000, clay_copula)
clay_copula <- claytonCopula(param = 0.5)
sample_copula2 <- rCopula(1000, clay_copula)
clay_copula <- claytonCopula(param = 1)
sample_copula3 <- rCopula(1000, clay_copula)
clay_copula <- claytonCopula(param = 2)
sample_copula4 <- rCopula(1000, clay_copula)

par(mfrow = c(2, 2))
plot(sample_copula1, main = "Independent")
plot(sample_copula2, main = expression(paste(theta, "= 0.5")))
plot(sample_copula3, main = expression(paste(theta, "= 1")))
plot(sample_copula4, main = expression(paste(theta, "= 2")))


## -----------------------------------------------------------------------------
A <- matrix(c(5, 4, 4, 5), ncol = 2)

results <- eigen(A, symmetric = TRUE,
                 only.values = TRUE)

c("GV" = prod(results$values), 
  "TV" = sum(results$values))

# Compare this with the following
B <- matrix(c(5, -4, -4, 5), ncol = 2)

# GV(A) = 9; TV(A) = 10
c("GV" = det(B), 
  "TV" = sum(diag(B)))

