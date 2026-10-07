
## Validation tests for daily rainfall grids ##

## Regression Kriging (MLE) ##

## Modelling rainfall as a sum of a linear trend with geographic
## covariates + a smooth Gaussian process with a spatial dependence

## This is the same as regression kriging, except both components
## are fit as a joint model using Maximum Likelihood Expectation


# %%

library(dplyr)

source("interpolation/MLE_kriging.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## Regression Formula
f <- as.formula(
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            # "dist2c", "exp25k",
            sep = " + "
        )
    )
)

# %%

rain_data_OK <- rain_data


for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        MLE_kriging(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],

            cov_function = "Exp",

            formula = f,

            init_pars = c(0.001, 0.0001, 150),
            lower = c(1e-8, 0, 10),
            upper = c(0.05, 0.01, 500)
        )$pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(paste("Date", date_index))
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "RK Exp MLE",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "RK_Exp_MLE_trend.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Check grids and computation time

times <- rep(NA, nrow(dates))

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        MLE_kriging(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            cov_function = "Exp",

            formula = f,

            init_pars = c(0.001, 0.0001, 150),
            lower = c(1e-8, 0, 10),
            upper = c(0.05, 0.01, 500)
        )$pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/RK_Exp_MLE_trend/RK_Exp_MLE_trend_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index] <- et - st

    print(date_index)
}

print("Mean Times")
print(mean(times, na.rm = TRUE))
