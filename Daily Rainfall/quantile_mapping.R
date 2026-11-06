
## Validation tests for daily rainfall grids ##

## Improving grids with quantile mapping ##

## As a final gridding step, I want to try quantile mapping ##
## This essentially means matching the distribution of the  ##
## gridded values with the distirbution of station values   ##

## I'm going to do it in a parametric way, using a fitted   ##
## Bernoulli-Gamma distribution to the station data (this   ##
## is well suited for Daily Rainfall)                       ##

## We'll apply it to one of my best gridding methods,       ##
## likelihood-based kriging                                 ##

## I'll keep this code around, but QM hasn't been very      ##
## effective thus far

# %%

## First do the basic setup
library(dplyr)
source("interpolation/IDW.R")
source("interpolation/MLE_kriging.R")
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

f <- y ~ east + north + points5

# %%

## Using the QMAP package for useful distributions
library(qmap)

## Define a quantile mapping function
quantile_map <- function(x1, x2) {

    if (sum(x2 > 0, na.rm = TRUE) < 5) {
        return(rep(0, length(x2)))
    }

    dists <- ifelse(c(any(x1 == 0), any(x2 == 0)), "berngamma", "gamma")

    pars <- list(
        gamma = c(
            "shape.o" = "shape",
            "rate.o" = "rate"
        ),
        berngamma = c(
            "prob.o" = "prob",
            "scale.o" = "scale",
            "shape.o" = "shape"
        )
    )

    x1_pars <- fitQmapDIST(x1, x1, distr = dists[1])$par
    x2_pars <- fitQmapDIST(x2, x2, distr = dists[2])$par

    getpars <- function(x, d) {
        p <- as.list(x[1, names(pars[[d]])])
        names(p) <- pars[[d]]
        p
    }

    funs <- list(gamma = list(pgamma, qgamma),
                 berngamma = list(pberngamma, qberngamma))

    p <- do.call(
        funs[[dists[2]]][[1]],
        c(list(q = x2), getpars(x2_pars, dists[2]))
    )

    do.call(
        funs[[dists[1]]][[2]],
        c(list(p = p), getpars(x1_pars, dists[1]))
    )
}

emp_quantile_map <- function(x1, x2) {

    p <- ecdf(x2)(x2)

    quantile(
        x1,
        probs = p,
        na.rm = TRUE,
        names = FALSE
    )
}


# %%

## Function that takes y_hat (predicted rain) and adds it to test set
update_daily_predictions <- function(
    rain_data,
    y_hat,
    dates,
    date_index,
    QM = FALSE
) {
    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    if (QM) {

        y_hat <-
            quantile_map(
                daily_rain_data$train$y,
                y_hat
            )
    }

    daily_rain_data$test <- daily_rain_data$test |>
        mutate(predicted_rain = expm1(y_hat) * LTA)

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    predicted_rain <-
        daily_rain_data$test$predicted_rain



    rain_data$test[
        rain_data$test$day == current_day &
            rain_data$test$month == current_month &
            rain_data$test$year == current_year,

        "predicted_rain"
    ] <-
        predicted_rain

    rain_data
}

# %%

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        IDW(
            daily_rain_data$train$y,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],

            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index,
        QM = TRUE
    )

    print(paste("Date", date_index))
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "RK Exp MLE QM",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)
