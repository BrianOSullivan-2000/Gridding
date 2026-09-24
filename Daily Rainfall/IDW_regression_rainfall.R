
## Validation tests for daily rainfall grids ##
## Inverse Distance Weighting with regression ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation, interpolate residuals using IDW ##

# %%

library(dplyr)
source("interpolation/regression_IDW.R")

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
            "dist2c", "exp25k",
            # "n5", "e5", "s5", "w5",
            # "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

# %%

## Stepwise Regression

source("regression/MLR.R")

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = MLR,
            formula = f,
            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW MLR",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_MLR.csv",
#     row.names = FALSE
# )


# %%

# Elastic-net Regularization

source("regression/elastic-net.R")

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = Elastic_Net,
            formula = f,
            idp = 2,
            nmax = 15,
            alpha = 0.9
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW ENET",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_ENET.csv",
#     row.names = FALSE
# )


# %%

# Generalized Least Squares

source("regression/GLS.R")

for (date_index in seq_len(nrow(dates))) {
    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = GLS,
            formula = f,
            flex_fit = TRUE, flex_vgm = TRUE,
            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW GLS",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_GLS.csv",
#     row.names = FALSE
# )

# %%

# Random Forests

source("regression/random_forest.R")

for (date_index in seq_len(nrow(dates))) {
    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],

            regression_method = Random_Forest,
            formula = f,

            mtry = 3,
            min.node.size = 1,

            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW RF",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_RF.csv",
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
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],

            regression_method = XGBOOST,
            formula = f,

            nrounds = 50,
            max_depth = 3,
            learning_rate = 0.01,
            min_child_weight = 1,
            subsample = 0.5,
            colsample_bytree = 0.5,

            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW XGB",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_XGB.csv",
#     row.names = FALSE
# )

# %%

## Plot grids and check computation time

times <- matrix(nrow = 200, ncol = 5)
colnames(times) <- c("Stepwise", "ENET", "GLS", "RF", "XGB")

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    ## MLR (Stepwise)

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_IDW(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            regression_method = MLR,
            formula = f,
            idp = 2,
            nmax = 15
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/IDW_MLR/IDW_MLR_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 1] <- et - st




    ## Elastic-net

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_IDW(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            regression_method = Elastic_Net,
            formula = f,
            idp = 2,
            nmax = 15,
            alpha = 0.9
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/IDW_ENET/IDW_ENET_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 2] <- et - st




    ## Generalised Least Squares

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_IDW(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            regression_method = GLS,
            formula = f,
            flex_fit = TRUE, flex_vgm = TRUE,
            idp = 2,
            nmax = 15
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/IDW_GLS/IDW_GLS_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 3] <- et - st




    ## Random Forest

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_IDW(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = Random_Forest,
            formula = f,

            mtry = 3,
            min.node.size = 1,
            cv_method = "oob",

            idp = 2,
            nmax = 15
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/IDW_RF/IDW_RF_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 4] <- et - st




    ## XGBoost

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_IDW(
            daily_rain_data,
            grid_geodata |> mutate(y = NA),
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = XGBOOST,
            formula = f,

            nrounds = 50,
            max_depth = 3,
            learning_rate = 0.01,
            min_child_weight = 1,
            subsample = 0.5,
            colsample_bytree = 0.5,

            idp = 2,
            nmax = 15
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/IDW_XGB/IDW_XGB_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 5] <- et - st

    print(date_index)
}

print(paste(
    "Mean Time",
    mean(times)
))
