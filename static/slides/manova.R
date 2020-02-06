## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache=FALSE)
options(knitr.kable.NA = '-')


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
n <- nrow(Y); p <- ncol(Y); g <- 2

W <- (nrow(Y_low) - 1)*cov(Y_low) + 
  (nrow(Y_high) - 1)*cov(Y_high)
B <- (n-1)*cov(Y) - W
(Lambda <- det(W)/det(W+B))


## -----------------------------------------------------------------------------
transf_lambda <- -(n - 1 - 0.5*(p + g))*log(Lambda)
transf_lambda > qchisq(0.95, p*(g-1))
# Or if you want a p-value
pchisq(transf_lambda, p*(g-1), lower.tail = FALSE)


## -----------------------------------------------------------------------------
# R has a function for MANOVA
# But first, create factor variable
rate <- gl(g, 10, labels = c("Low", "High"))

fit <- manova(Y ~ rate)
summary_tbl <- broom::tidy(fit, test = "Wilks")
# Or you can use the summary function


## -----------------------------------------------------------------------------
knitr::kable(summary_tbl, digits = 3)


## ----message = FALSE----------------------------------------------------------
# Check residuals for evidence of normality
library(tidyverse)
resids <- residuals(fit)

data_plot <- gather(as.data.frame(resids),
                    variable, residual)

ggplot(data_plot, aes(sample = residual)) + 
  stat_qq() + stat_qq_line() +
  facet_grid(. ~ variable) + 
  theme_minimal()


## -----------------------------------------------------------------------------
# Next: Chi-squared plot
Sn <- cov(resids)
dists <- mahalanobis(resids, colMeans(resids), Sn)
df <- mean(dists)

qqplot(qchisq(ppoints(dists), df = df),
       dists, xlab = "Theoretical Quantiles")
qqline(dists, distribution = function(p) {
  qchisq(p, df = df)
  })


## -----------------------------------------------------------------------------
knitr::kable(broom::tidy(fit), digits = 3)


## -----------------------------------------------------------------------------
knitr::kable(broom::tidy(fit, test = "Hotelling-Lawley"), 
             digits = 3)

