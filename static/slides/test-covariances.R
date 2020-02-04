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
df <- choose(p + 1, 2)
c(lrt, qchisq(0.95, df))


## -----------------------------------------------------------------------------
lrt <- -2*0.5*n*(log(det(V)) - p*log(mean(diag(V))))
df <- choose(p + 1, 2) - 1

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


## -----------------------------------------------------------------------------
url <- paste0("https://maxturgeon.ca/w20-stat7200/",
              "blue_data.csv")
blue_data <- read.csv(url)
names(blue_data)
dim(blue_data)


## -----------------------------------------------------------------------------
# Let's test for independence between 
# all four variables
n <- nrow(blue_data)
p <- ncol(blue_data)

V <- (n-1)*cov(blue_data)
lrt <- -2*(log(det(V)) - sum(log(diag(V))))


## -----------------------------------------------------------------------------
df <- choose(p + 1, 2) - p
c(lrt, qchisq(0.95, df))
lrt > qchisq(0.95, df)


## -----------------------------------------------------------------------------
## Example on producing plastic film 
## from Krzanowski (1998, p. 381)
tear <- c(6.5, 6.2, 5.8, 6.5, 6.5, 6.9, 7.2, 
          6.9, 6.1, 6.3, 6.7, 6.6, 7.2, 7.1, 
          6.8, 7.1, 7.0, 7.2, 7.5, 7.6)
gloss <- c(9.5, 9.9, 9.6, 9.6, 9.2, 9.1, 10.0, 
           9.9, 9.5, 9.4, 9.1, 9.3, 8.3, 8.4, 
           8.5, 9.2, 8.8, 9.7, 10.1, 9.2)
opacity <- c(4.4, 6.4, 3.0, 4.1, 0.8, 5.7, 2.0, 
             3.9, 1.9, 5.7, 2.8, 4.1, 3.8, 1.6, 
             3.4, 8.4, 5.2, 6.9, 2.7, 1.9)


## -----------------------------------------------------------------------------
Y <- cbind(tear, gloss, opacity)
Y_low <- Y[1:10,]
Y_high <- Y[11:20,]
n <- nrow(Y); p <- ncol(Y); K <- 2
n1 <- n2 <- nrow(Y_low)


## -----------------------------------------------------------------------------
Sig_low <- (n1 - 1)*cov(Y_low)/n1
Sig_high <- (n2 - 1)*cov(Y_high)/n2
Sig_pool <- (n1*Sig_low + n2*Sig_high)/n

c("pool" = log(det(Sig_pool)),
  "low" = log(det(Sig_low)),
  "high" = log(det(Sig_high)))


## -----------------------------------------------------------------------------
lrt <- n*log(det(Sig_pool)) - 
  n1*log(det(Sig_low)) - 
  n2*log(det(Sig_high))
df <- (K - 1)*choose(p + 1, 2)
c(lrt, qchisq(0.95, df))


## -----------------------------------------------------------------------------
S_low <- cov(Y_low)
S_high <- cov(Y_high)
S_pool <- ((n1 - 1)*S_low + (n2 - 1)*S_high)/(n - K)

lrt2 <- (n - K)*log(det(S_pool)) - 
  (n1 - 1)*log(det(S_low)) - 
  (n2 - 1)*log(det(S_high))

c(lrt, lrt2, qchisq(0.95, df))


## -----------------------------------------------------------------------------
u <- (2*p^2 + 3*p - 1)/(6*(p + 1)*(K - 1))
u <- u * ((n1 - 1)^{-1} + (n2 - 1)^{-1} - (n - K)^{-1})
lrt3 <- lrt2*(1 - u)

c(lrt, lrt2, lrt3, qchisq(0.95, df))


## ----message = FALSE----------------------------------------------------------
# You can also visualize the covariances----
library(heplots)
rate <- gl(K, 10, labels = c("Low", "High"))
boxm_res <- boxM(Y, rate)

# You can plot the log generalized variances
# The plot function adds 95% CI
plot(boxm_res)


## -----------------------------------------------------------------------------
# Finally you can also plot the ellipses
# as a way to compare the covariances
covEllipses(Y, rate, center = TRUE, 
            label.pos = 'bottom')


## -----------------------------------------------------------------------------
# Or all pairwise comparisons together
covEllipses(Y, rate, center = TRUE, 
            label.pos = 'bottom', 
            variables = 1:3)


## -----------------------------------------------------------------------------
set.seed(7200)

# Simulation parameters
n <- 10
p <- 2
B <- 1000


## -----------------------------------------------------------------------------
# Generate data
lrt_dist <- replicate(B, {
  Y <- matrix(rnorm(n*p), ncol = p)
  V <- crossprod(Y)
  # log Lambda
  0.5*n*(log(det(V)) - p*log(mean(diag(V))))
})


## -----------------------------------------------------------------------------
# General asymptotic result
df <- choose(p + 1, 2)
general_chisq <- rchisq(B, df = df)


## -----------------------------------------------------------------------------
# Bartlett's correction
df <- choose(p + 1, 2) - 1
const <- (6*p*(n-1) - (2*p^2 + p + 2))/(6*p*n)
bartlett_chisq <- rchisq(B, df = df)/const


## -----------------------------------------------------------------------------
# Plot empirical CDFs
plot(ecdf(-2*lrt_dist), main = "-2 log Lambda")
lines(ecdf(general_chisq), col = 'blue')
lines(ecdf(bartlett_chisq), col = 'red')
legend('bottomright', 
       legend = c("-2log Lambda", "General approx.", 
                  "Bartlett"),
       lty = 1, col = c('black', 'blue', 'red'))

