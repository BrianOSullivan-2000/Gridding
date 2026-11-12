
## Validation tests for daily rainfall grids ##

## Ordinary Kriging ##
## Assumes data follows a Gaussian process with a spatial dependence ##
## Data is predicted according to distance from nearby points ##

# %%

library(dplyr)
source("interpolation/ordinary_kriging.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_provisional.RData")

# %%

## Set up grid search for hyperparameters
cutoffs <- c(150000, 200000, 250000, 300000)
widths <- c(10000, 15000, 20000, 25000)

hyperparameters <- expand.grid(cutoff = cutoffs, width = widths)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    cutoff <- hyperparameters[hyperparameter_index, ]$cutoff
    width <- hyperparameters[hyperparameter_index, ]$width

    for (date_index in seq_len(nrow(dates))) {

        daily_rain_data <- get_daily_rain_data(
            rain_data, dates, date_index
        )

        y_hat <-
            ordinary_kriging(
                daily_rain_data$train$y,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],

                cutoff = cutoff,
                width = width,
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
            metrics_table,
            paste0("OK cutoff-", cutoff, " width-", width),
            rain_data$test$rain,
            rain_data$test$predicted_rain
        )

    # print(metrics_table)
}

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_provisional/",
#         "OK_Exp_config_log.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Also going to check against nearest neighbours
nmaxs <- c(8, 12, 15, 20, 50, Inf)
cutoffs <- c(200000, 250000, 300000)

hyperparameters <- expand.grid(cutoff = cutoffs, nmax = nmaxs)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    cutoff <- hyperparameters[hyperparameter_index, ]$cutoff
    nmax <- hyperparameters[hyperparameter_index, ]$nmax

    for (date_index in seq_len(nrow(dates))) {

        daily_rain_data <- get_daily_rain_data(
            rain_data, dates, date_index
        )

        y_hat <-
            ordinary_kriging(
                daily_rain_data$train$y,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],

                cutoff = cutoff,
                nmax = nmax,

                width = 15000,
                flex_fit = TRUE,
                vgm_model = "Exp",
                kappa = 0.5,
                debug.level = 0
            )$var1.pred

        rain_data <- update_daily_predictions(
            rain_data, y_hat, dates, date_index
        )
    }

    metrics_table <-
        collect_metrics(
            metrics_table,
            paste0("OK cutoff-", cutoff, " nmax-", nmax),
            rain_data$test$rain,
            rain_data$test$predicted_rain
        )

    # print(metrics_table)
}

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_provisional/",
#         "OK_Exp_config_nmax_log.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Check grids and computation time

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        ordinary_kriging(
            daily_rain_data$y,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            cutoff = 250000,
            width = 15000,
            nmax = 15,

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
            "Figures/Daily_Rainfall_provisional/OK_Exp/OK_Exp_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    print(date_index)
}
