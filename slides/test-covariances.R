## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=FALSE)


## ---- message = FALSE---------------------------------------------------------
library(tidyverse)
# Winnipeg avg temperature
url <- paste0("https://maxturgeon.ca/w20-stat7200/",
              "winnipeg_temp.csv")
dataset <- read.csv(url)
dataset[1:3,1:3]


## -----------------------------------------------------------------------------
n <- nrow(dataset)
p <- ncol(dataset)

V <- (n - 1)*cov(dataset)


## -----------------------------------------------------------------------------
# Diag = 14^2
# Corr = 0.8
Sigma0 <- diag(0.8, nrow = p)
diag(Sigma0) <- 1
Sigma0 <- 14^2*Sigma0
Sigma0_invXV <- solve(Sigma0, V)


## -----------------------------------------------------------------------------
lrt <- 0.5*n*p*(1 - log(n))
lrt <- lrt + 0.5*n*log(det(Sigma0_invXV))
lrt <- lrt - 0.5*sum(diag(Sigma0_invXV))
lrt <- -2*lrt


## -----------------------------------------------------------------------------
df <- choose(p, 2)
c(lrt, qchisq(0.95, df))


## -----------------------------------------------------------------------------
lrt <- -2*0.5*n*(log(det(V)) - p*log(mean(diag(V))))
df <- choose(p, 2) - 1

c(lrt, qchisq(0.95, df))


## -----------------------------------------------------------------------------
B <- 1000
df1 <- 0.5*(n - seq_len(p-1) - 1)
df2 <- seq_len(p-1)*(0.5 + 1/p)

# Critical values
dist <- replicate(B, {
  prod(rbeta(p-1, df1, df2))
  })


## -----------------------------------------------------------------------------
# Test statistic
decomp <- eigen(V, symmetric = TRUE, only.values = TRUE)
ar_mean <- mean(decomp$values)
geo_mean <- exp(mean(log(decomp$values)))

lrt_mod <- (geo_mean/ar_mean)^p

c(lrt_mod, quantile(dist, 0.95))

