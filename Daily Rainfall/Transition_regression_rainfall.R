
## Validation tests for daily rainfall grids ##
## Transition Regression by Reto Stauffer    ##

## Regression-based approach that assumes daily rainfall follows a
## transition model, and fits it using a smooth relationship with
## geographic covariates

# %%

library(dplyr)

source("interpolation/regression_grid.R")
source("regression/transitreg.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")


# %%

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    if (length(unique(daily_rain_data$train$y)) < 10) {
        k_theta = 2
    } else {
        k_theta = 5
    }

    ## Define formula (include censoring if zeros present)
    if (sum(daily_rain_data$train$y >  0.001) < 10) {
        censored <- "left"

        transitreg_formula <- as.formula(
            substitute(
                y ~
                    theta0 +
                        s(east, k = 10) +
                        s(north, k = 10) +
                        ti(east, north, k = 20) +
                        s(points5, k = 10),
                list(K = k_theta)
            )
        )
    } else if (any(daily_rain_data$train$rain) == 0) {
        censored <- "left"

        transitreg_formula <- as.formula(
            substitute(
                y ~
                    theta0 +
                        s(theta, k = K) +
                        s(east, k = 10) +
                        s(north, k = 10) +
                        ti(east, north, k = 20) +
                        s(points5, k = 10),
                list(K = k_theta)
            )
        )

    } else {
        censored <- "uncensored"

        transitreg_formula <- as.formula(
            substitute(
                y ~
                    s(theta, k = K) +
                        s(east, k = 10) +
                        s(north, k = 10) +
                        ti(east, north, k = 20) +
                        s(points5, k = 10),
                list(K = k_theta)
            )
        )
    }

    breaks <- seq(
        0,
        ceiling(max(daily_rain_data$train$y) / 0.01) * 0.01,
        by = 0.01
    )
    if (length(breaks) < 10) {
        breaks <-
            seq(
                0,
                ceiling(max(daily_rain_data$train$y) / 0.01) * 0.01,
                length.out = 10
            )
    }

    y_hat <-
        regression_grid(
            daily_rain_data$train,
            daily_rain_data$test,
            regression_method = TransitReg,
            formula = transitreg_formula,
            breaks = breaks,
            censored = censored
        )$pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "Transitreg",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/Transitreg.csv",
#     row.names = FALSE
# )

# %%


## Plot grids and check computation time

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

# for (date_index in seq_len(nrow(dates))) {


## Crashed for plot 65 (need to run again) ##
x <- seq_len(nrow(dates))
for (date_index in x[x > 64]) {




    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index,
        experiment_type = "all_data"
    )

    if (length(unique(daily_rain_data$y)) < 10) {
        k_theta = 2
    } else {
        k_theta = 5
    }

    ## Define formula (include censoring if zeros present)
    if (sum(daily_rain_data$y >  0.001) < 10) {
        censored <- "left"

        transitreg_formula <- as.formula(
            substitute(
                y ~
                    theta0 +
                        s(east, k = 10) +
                        s(north, k = 10) +
                        ti(east, north, k = 20) +
                        s(points5, k = 10),
                list(K = k_theta)
            )
        )
    } else if (any(daily_rain_data$rain) == 0) {
        censored <- "left"

        transitreg_formula <- as.formula(
            substitute(
                y ~
                    theta0 +
                        s(theta, k = K) +
                        s(east, k = 10) +
                        s(north, k = 10) +
                        ti(east, north, k = 20) +
                        s(points5, k = 10),
                list(K = k_theta)
            )
        )

    } else {
        censored <- "uncensored"

        transitreg_formula <- as.formula(
            substitute(
                y ~
                    s(theta, k = K) +
                        s(east, k = 10) +
                        s(north, k = 10) +
                        ti(east, north, k = 20) +
                        s(points5, k = 10),
                list(K = k_theta)
            )
        )
    }

    breaks <- seq(
        0,
        ceiling(max(daily_rain_data$y) / 0.01) * 0.01,
        by = 0.01
    )
    if (length(breaks) < 10) {
        breaks <-
            seq(
                0,
                ceiling(max(daily_rain_data$y) / 0.01) * 0.01,
                length.out = 10
            )
    }

    y_hat <-
        regression_grid(
            daily_rain_data,
            grid_geodata,
            regression_method = TransitReg,
            formula = transitreg_formula,
            breaks = breaks,
            censored = censored
        )$pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/Transitreg/Transitreg_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)

    print(date_index)
}

print(paste(
    "Mean Time",
    mean(times)
))
