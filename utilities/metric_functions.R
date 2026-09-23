
## Source file for several Metric functions ##

## Metrics package has Root Mean Squared Error (among others)
library(Metrics)

## Coefficient of Determination
R2 <- function(pred, obs) {
    1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

## Mean Error
me <- function(pred, obs) {
    mean(pred - obs)
}

## Jensen-Shannon Divergence (compares empirical distributions)
JSD <- function(y, y_hat) {
    p <- y / sum(y)
    q <- y_hat / sum(y_hat)

    m <- 0.5 * (p + q)
    kl_p_m <- sum(ifelse(p > 0, p * log2(p / m), 0))
    kl_q_m <- sum(ifelse(q > 0, q * log2(q / m), 0))

    JSD <- 0.5 * kl_p_m + 0.5 * kl_q_m
    JSD
}

## Ten largest errors
ME10 <- function(y_hat, y) {
    errs <- round(y_hat - y, 2)
    max_errs <- paste(sort(errs, decreasing = TRUE)[1:10], collapse = ", ")
    max_errs
}

## Check bias for Wet Days (Daily Rainfall)
Wet_Day_Bias <- function(y_hat, y) {
    (sum(y_hat > 0.2) - sum(y > 0.2)) * 100 / sum(y > 0.2)
}