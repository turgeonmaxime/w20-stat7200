## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## -----------------------------------------------------------------------------
Sigma_Y <- matrix(c(1, 0.4, 0.4, 1), ncol = 2)
Sigma_X <- matrix(c(1, 0.2, 0.2, 1), ncol = 2)
Sigma_YX <- matrix(c(0.5, 0.3, 0.6, 0.4), ncol = 2)
Sigma_XY <- t(Sigma_YX)

rbind(cbind(Sigma_Y, Sigma_YX),
      cbind(Sigma_XY, Sigma_X))

## ---- message = FALSE---------------------------------------------------------
library(expm)
sqrt_Y <- sqrtm(Sigma_Y)
sqrt_X <- sqrtm(Sigma_X)
M1 <- solve(sqrt_Y) %*% Sigma_YX %*% solve(Sigma_X)%*% 
  Sigma_XY %*% solve(sqrt_Y)

(decomp1 <- eigen(M1))


## -----------------------------------------------------------------------------
decomp1$vectors[,1] %*% solve(sqrt_Y)


## -----------------------------------------------------------------------------
M2 <- solve(sqrt_X) %*% Sigma_XY %*% solve(Sigma_Y)%*% 
  Sigma_YX %*% solve(sqrt_X)

decomp2 <- eigen(M2)
decomp2$vectors[,1] %*% solve(sqrt_X)


## -----------------------------------------------------------------------------
sqrt(decomp1$values)


## -----------------------------------------------------------------------------
# Let's generate data
library(mvtnorm)
Sigma <- rbind(cbind(Sigma_Y, Sigma_YX),
               cbind(Sigma_XY, Sigma_X))

YX <- rmvnorm(100, sigma = Sigma)
Y <- YX[,1:2]
X <- YX[,3:4]

decomp <- stats::cancor(x = X, y = Y)


## -----------------------------------------------------------------------------
U <- Y %*% decomp$ycoef
V <- X %*% decomp$xcoef

diag(cor(U, V))
decomp$cor


## ---- message=FALSE-----------------------------------------------------------
library(tidyverse)
library(dslabs)

str(olive)


## ---- message=FALSE-----------------------------------------------------------
# X contains the type of acids
X <- select(olive, -area, -region) %>% 
  as.matrix

# Y contains the information about regions
count(olive, region)
Y <- select(olive, region) %>% 
  model.matrix(~ region - 1, data = .)


## ---- message=FALSE-----------------------------------------------------------
# We get three dummy variables
head(unname(Y))

decomp <- cancor(X, Y)

V <- X %*% decomp$xcoef


## ---- message=FALSE-----------------------------------------------------------
data.frame(
  V1 = V[,1],
  V2 = V[,2],
  region = olive$region
) %>% 
  ggplot(aes(V1, V2, colour = region)) +
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'top')


## ---- echo = -1---------------------------------------------------------------
old_opts <- options(digits = 2)
# Olive data--Standardize
X_sc <- scale(X)
Y_sc <- scale(Y)
decomp_sc <- cancor(X_sc, Y_sc)

# Extract Canonical variates
V_sc <- X_sc %*% decomp_sc$xcoef
colnames(V_sc) <- paste0("CC", seq_len(ncol(V_sc)))


## ---- echo = -1---------------------------------------------------------------
old_opts <- options(digits = 2)
(prop_X <- rowMeans(cor(V_sc, X_sc)^2))

cumsum(prop_X)


## ---- echo = -1---------------------------------------------------------------
old_opts <- options(digits = 2)
# But since we are dealing with correlations
# We get the same with unstandardized variables
decomp <- cancor(X, Y)
V <- X %*% decomp$xcoef
colnames(V) <- paste0("CC", seq_len(ncol(V)))

(prop_X <- rowMeans(cor(V, X)^2))

cumsum(prop_X)


## ---- echo = -1---------------------------------------------------------------
options(digits = old_opts$digits)
# Let's go back to the olive data
decomp <- cancor(X, Y)
V <- X %*% decomp$xcoef
colnames(V) <- paste0("CC", seq_len(8))

library(lattice)
levelplot(cor(X, V[,1:2]), 
          at = seq(-1, 1, by = 0.1),
          xlab = "", ylab = "")


## -----------------------------------------------------------------------------
levelplot(cor(Y, V[,1:2]), 
          at = seq(-1, 1, by = 0.1),
          xlab = "", ylab = "")


## -----------------------------------------------------------------------------
# Canonical correlations
decomp$cor

# Maximum value in correlation matrix
max(abs(cor(Y, X)))

# Multiple correlation coefficients
sqrt(summary(lm(V[,1] ~ Y))$r.squared)
sqrt(summary(lm(V[,2] ~ Y))$r.squared)


## ---- message = FALSE, echo = FALSE-------------------------------------------
Sigma_Y <- matrix(c(1, 0.4, 0.4, 1), ncol = 2)
Sigma_X <- matrix(c(1, 0.2, 0.2, 1), ncol = 2)
Sigma_YX <- matrix(c(0.5, 0.3, 0.6, 0.4), ncol = 2)
Sigma_XY <- t(Sigma_YX)

library(expm)
sqrt_Y <- sqrtm(Sigma_Y)
sqrt_X <- sqrtm(Sigma_X)
M1 <- solve(sqrt_Y) %*% Sigma_YX %*% solve(Sigma_X)%*% 
  Sigma_XY %*% solve(sqrt_Y)

decomp1 <- eigen(M1, symmetric = TRUE)
decompY <- eigen(Sigma_Y, symmetric = TRUE)


## ----echo = FALSE-------------------------------------------------------------
plot(x = 0, y = 0, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "Y1", ylab = "Y2", type = 'n', asp = 1)
arrows(x0 = 0, y0 = 0,
       x1 = 1, y1 = 0, col = 'red')
arrows(x0 = 0, y0 = 0,
       x1 = 0, y1 = 1, col = 'red')


## ----echo = FALSE-------------------------------------------------------------
Ytrans <- decompY$vectors
plot(x = 0, y = 0, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "Y1", ylab = "Y2", type = 'n', asp = 1)
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,1], y1 = Ytrans[2,1], 
       col = 'red')
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,2], y1 = Ytrans[2,2], 
       col = 'red')


## ----echo = FALSE-------------------------------------------------------------
Ytrans <- decompY$vectors %*% diag(1/sqrt(decompY$values))
plot(x = 0, y = 0, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "Y1", ylab = "Y2", type = 'n', asp = 1)
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,1], y1 = Ytrans[2,1], 
       col = 'red')
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,2], y1 = Ytrans[2,2], 
       col = 'red')


## ----echo = FALSE-------------------------------------------------------------
Ytrans <- decompY$vectors %*% (decompY$vectors %*% diag(1/sqrt(decompY$values)))
plot(x = 0, y = 0, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "Y1", ylab = "Y2", type = 'n', asp = 1)
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,1], y1 = Ytrans[2,1], 
       col = 'red')
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,2], y1 = Ytrans[2,2], 
       col = 'red')


## ----echo = FALSE-------------------------------------------------------------
Ytrans <- decomp1$vectors %*% solve(sqrt_Y)
plot(x = 0, y = 0, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "Y1", ylab = "Y2", type = 'n', asp = 1)
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,1], y1 = Ytrans[2,1], 
       col = 'red')
arrows(x0 = 0, y0 = 0,
       x1 = Ytrans[1,2], y1 = Ytrans[2,2], 
       col = 'red')


## -----------------------------------------------------------------------------
library(vegan)

data(varespec)
data(varechem)

# There are too many variables in varespec
# Let's pick first 10
Y <- select(varespec, Callvulg:Diphcomp) %>% 
  as.matrix


## -----------------------------------------------------------------------------
# The help page in `vegan` suggests a better 
# chemical model
X <- model.matrix( ~ Al + P*(K + Baresoil) - 1,
                  data = varechem)
colnames(X)[1:4]
colnames(X)[5:6]


## -----------------------------------------------------------------------------
decomp <- cancor(x = X, y = Y)

n <- nrow(X)
(LRT <- -n*log(prod(1 - decomp$cor^2)))

p <- min(ncol(X), ncol(Y))
q <- max(ncol(X), ncol(Y))
LRT > qchisq(0.95, df = p*q)

LRT_bart <- -(n - 1 - 0.5*(p + q + 1)) *
  log(prod(1 - decomp$cor^2))

c("Large Sample" = LRT,
  "Bartlett" = LRT_bart)

LRT_bart > qchisq(0.95, df = p*q)


## ----echo = -1----------------------------------------------------------------
old_opts <- options(digits = 2)
# We can get the truncated LRTs in one go
(log_ccs <- rev(log(cumprod(1 - rev(decomp$cor)^2))))

(LRTs <- -(n - 1 - 0.5*(p + q + 1)) * log_ccs)

k_seq <- seq(0, p - 1)
LRTs > qchisq(0.95,
              df = (p - k_seq)*(q - k_seq))
# We only reject the first null hypothesis 
# of independence


## ----echo = FALSE-------------------------------------------------------------
data.frame(
  k = k_seq,
  LRT = LRTs,
  CC = qchisq(0.95,
              df = (p - k_seq)*(q - k_seq))
) %>% 
  gather(Type, Value, LRT, CC) %>% 
  ggplot(aes(k, Value, colour = Type)) + 
  geom_point() +
  geom_line() +
  theme_minimal() +
  ylab("") +
  scale_colour_discrete(name = "",
                        breaks = c("LRT", "CC"),
                        labels = c("LRT", "Critical value"))


## ---- message = FALSE---------------------------------------------------------
# Recall the plastic film data
library(heplots)

fit <- lm(cbind(tear, gloss, opacity) ~ rate + additive,
          data = Plastic)
coef(fit)


## -----------------------------------------------------------------------------
Y <- Plastic %>% 
  select(tear, gloss, opacity) %>% 
  as.matrix
X <- model.matrix(~ rate + additive, data = Plastic)

# We get the same as OLS
(beta_ols <- solve(crossprod(X), crossprod(X, Y)))


## -----------------------------------------------------------------------------
# Reduced-Rank regression
M <- crossprod(Y, X) %*% beta_ols
decomp <- eigen(M)

# Take rank = 1
W <- decomp$vectors[,1, drop=FALSE]
rownames(W) <- colnames(Y)
(beta_rrr <- beta_ols %*% tcrossprod(W))

# Note that rank 1 means rows are colinear
beta_rrr[1,]/beta_rrr[2,]


## -----------------------------------------------------------------------------
# Let's create a function
redrank <- function(Y, X, rank = 1) {
  beta_ols <- solve(crossprod(X), crossprod(X, Y))
  M <- crossprod(Y, X) %*% beta_ols
  decomp <- eigen(M)
  W <- decomp$vectors[,seq_len(rank),drop=FALSE]
  rownames(W) <- colnames(Y)
  return(beta_ols %*% tcrossprod(W))
}


## -----------------------------------------------------------------------------
all.equal(beta_rrr, redrank(Y, X))


## -----------------------------------------------------------------------------
# First the log likelihoods
loglik <- sapply(c(1, 2, 3), function(k) {
  beta_rrr <- redrank(Y, X, k)
  resids <- Y - X %*% beta_rrr
  n*log(det(crossprod(resids)/nrow(Y)))
})


## -----------------------------------------------------------------------------
# With naive degrees of freedom
2*seq_len(3)*(ncol(X) + ncol(Y) - 
                seq_len(3)) + loglik


## -----------------------------------------------------------------------------
# With exact degrees of freedom
dfs <- sapply(seq_len(3), function(k) {
  total <- 0
  lambdas <- decomp$values[seq(k+1, ncol(Y))]
  for (ell in seq(1, k)) {
    total <- sum(lambdas/(decomp$values[ell] - lambdas)) + total
  }
  if (k == ncol(Y)) return(0) else return(2*total)
})


## -----------------------------------------------------------------------------
2*seq_len(3)*(ncol(X) + ncol(Y) - 
                seq_len(3)) + 2*dfs + loglik
# Both approaches select the full rank model


## -----------------------------------------------------------------------------
# Constrast this with rrpack::rrr
# Which uses a different AIC
rrpack::rrr(Y, X, ic.type = "AIC")


## -----------------------------------------------------------------------------
# Tobacco dataset
tobacco_y <- as.matrix(rrr::tobacco[,1:3])
tobacco_x <- as.matrix(rrr::tobacco[,4:9])

dim(tobacco_x)
dim(tobacco_y)


## -----------------------------------------------------------------------------
(rr_fit <- rrpack::rrr(tobacco_y, tobacco_x))


## -----------------------------------------------------------------------------
library(lattice)
coef <- rr_fit$coef
colnames(coef) <- colnames(tobacco_y)
levelplot(coef)

