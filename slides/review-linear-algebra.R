## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)


## -----------------------------------------------------------------------------
(A <- matrix(c(1, 2, 3, 2), ncol = 2))

eigen(A)$values


## -----------------------------------------------------------------------------
eigen(A)$vectors


## -----------------------------------------------------------------------------
A <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

# Eigenvalue method
result <- eigen(A)
Lambda <- diag(result$values)
P <- result$vectors
A_sqrt <- P %*% Lambda^0.5 %*% t(P)

all.equal(A, A_sqrt %*% A_sqrt) # CHECK

# Cholesky method
# It's upper triangular!
(L <- chol(A))

all.equal(A, t(L) %*% L) # CHECK


## ---- cache = FALSE-----------------------------------------------------------
set.seed(1234)
A <- matrix(rnorm(3 * 2), ncol = 2, nrow = 3)
result <- svd(A)
names(result)

result$d
result$u
result$v

D <- diag(result$d)
all.equal(A, result$u %*% D %*% t(result$v)) #CHECK


## ---- cache = FALSE-----------------------------------------------------------
# Note: crossprod(A) == t(A) %*% A
# tcrossprod(A) == A %*% t(A)
U <- eigen(tcrossprod(A))$vectors
V <- eigen(crossprod(A))$vectors

D <- matrix(0, nrow = 3, ncol = 2)
diag(D) <- result$d

all.equal(A, U %*% D %*% t(V)) # CHECK


## ---- cache = FALSE-----------------------------------------------------------
# What went wrong?
# Recall that eigenvectors are unique 
# only up to a sign!

# These elements should all be positive
diag(t(U) %*% A %*% V)

# Therefore we need to multiply the 
# corresponding columns of U or V 
# (but not both!) by -1
cols_flip <- which(diag(t(U) %*% A %*% V) < 0)
V[,cols_flip] <- -V[,cols_flip]

all.equal(A, U %*% D %*% t(V)) # CHECK

