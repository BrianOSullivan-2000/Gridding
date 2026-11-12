
## Validation tests for daily rainfall grids ##

## Regression Kriging ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation with kriging, assuminig regression residuals   ##
##  follow a Gaussian Process (with spatial dependence)                   ##


# %%

library(dplyr)

source("interpolation/regression_kriging.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_provisional.RData")

# %%

## Regression Formula
f_simple <- as.formula(
    paste(
        "y ~",
        paste(
            "east", "north",
            # "points5",
            sep = " + "
        )
    )
)
f_complex <- as.formula(
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            "dist2c", "exp25k",
            "n5", "e5", "s5", "w5",
            "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

f <- f_simple

# %%

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_kriging(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = MLR,
            formula = f_simple,

            cutoff = 250000,
            width = 15000,
            nmax = 12,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "RK Exp MLR",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_provisional/",
#         "RK_Exp.csv"
#     ),
#     row.names = FALSE
# )


# %%

# XGBoost

source("regression/xgboost.R")

for (date_index in seq_len(nrow(dates))) {
    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_kriging(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],

            regression_method = XGBOOST,
            formula = f_complex,

            nrounds = 300,
            max_depth = 6,
            learning_rate = 0.001,
            min_child_weight = 1,
            subsample = 0.2,
            colsample_bytree = 1,

            cutoff = 250000,
            width = 15000,
            nmax = 12,
            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0

        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )

    if (date_index %% 20 == 0) {
        print(date_index)
    }
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "RK XGB",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_provisional/",
#         "RK_Exp_XGB.csv"
#     ),
#     row.names = FALSE
# )


# %%

## Check grids and computation time

## Checking Exp with MLR (simple model)
## and Exp with XGBoost (complex model)

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    ## Exp MLR

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_kriging(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = MLR,
            formula = f_simple,

            cutoff = 250000,
            width = 15000,
            nmax = 12,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall_provisional/RK_Exp/RK_Exp_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )


    ## Exp XGBoost

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_kriging(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = XGBOOST,
            formula =f_complex,

            nrounds = 300,
            max_depth = 6,
            learning_rate = 0.001,
            min_child_weight = 1,
            subsample = 0.2,
            colsample_bytree = 1,

            cutoff = 250000,
            width = 15000,
            nmax = 12,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall_provisional/RK_Exp_XGB/RK_Exp_XGB_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    print(date_index)
}