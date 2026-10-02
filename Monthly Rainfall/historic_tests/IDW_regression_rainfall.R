
## Validation tests for monthly rainfall grids ##
## Inverse Distance Weighting with regression ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation, interpolate residuals using IDW ##

# %%

library(dplyr)

source("interpolation/regression_IDW.R")

## Load starting data
load("Data/Monthly_Rainfall/train_test_80_20_1883-1940.RData")

# %%

## Regression Formula
f <- as.formula(
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            # "dist2c",
            # "exp25k",
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

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            monthly_rain_data$train,
            monthly_rain_data$test,
            monthly_rain_data$train[c("east", "north")],
            monthly_rain_data$test[c("east", "north")],
            regression_method = MLR,
            formula = f,
            idp = 1.5,
            nmax = 15
        )$var1.pred

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW MLR",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Monthly_Rainfall/train_test_80_20_1883-1940/IDW_MLR.csv",
#     row.names = FALSE
# )

# %%

# Elastic-net Regularization

source("regression/elastic-net.R")

for (date_index in seq_len(nrow(dates))) {

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            monthly_rain_data$train,
            monthly_rain_data$test,
            monthly_rain_data$train[c("east", "north")],
            monthly_rain_data$test[c("east", "north")],
            regression_method = Elastic_Net,
            formula = f,
            idp = 1.5,
            nmax = 15,
            alpha = 0.1
        )$var1.pred

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW ENET",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Monthly_Rainfall/train_test_80_20_1883-1940/IDW_ENET.csv",
#     row.names = FALSE
# )


# %%

## Tried GLS but no luck - singular solutions
## I suspect it might be because overlapping stations?

## I don't expect IDW GLS to be especially good or anything
## So leaving it out for now

# %%

# Random Forests

source("regression/random_forest.R")

for (date_index in seq_len(nrow(dates))) {
    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            monthly_rain_data$train,
            monthly_rain_data$test,
            monthly_rain_data$train[c("east", "north")],
            monthly_rain_data$test[c("east", "north")],

            regression_method = Random_Forest,
            formula = f,

            mtry = 3,
            min.node.size = 1,

            idp = 1.5,
            nmax = 15
        )$var1.pred

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
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
#     "Results/Monthly_Rainfall/train_test_80_20_1883-1940/IDW_RF.csv",
#     row.names = FALSE
# )

# %%

# XGBoost

source("regression/xgboost.R")

for (date_index in seq_len(nrow(dates))) {
    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            monthly_rain_data$train,
            monthly_rain_data$test,
            monthly_rain_data$train[c("east", "north")],
            monthly_rain_data$test[c("east", "north")],

            regression_method = XGBOOST,
            formula = f,

            nrounds = 100,
            max_depth = 6,
            learning_rate = 0.001,
            min_child_weight = 1,
            subsample = 0.2,
            colsample_bytree = 1,

            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
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
#     "Results/Monthly_Rainfall/train_test_80_20_1883-1940/IDW_XGB.csv",
#     row.names = FALSE
# )

# %%

## Plot grids and check computation time

times <- matrix(nrow = nrow(dates), ncol = 2)
colnames(times) <- c("MLR", "XGB")

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    ## Generalised Least Squares

    st <- Sys.time()

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_IDW(
            monthly_rain_data,
            grid_geodata,
            monthly_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            regression_method = MLR,
            formula = f,
            idp = 1.5,
            nmax = 15
        )$var1.pred

    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    monthly_rain_plot(
        grid_geodata,
        monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/historic/IDW_MLR/IDW_MLR_",
            current_year, "_",
            sprintf("%02d", current_month),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 1] <- et - st




    ## XGBoost

    st <- Sys.time()

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_IDW(
            monthly_rain_data,
            grid_geodata |> mutate(y = NA),
            monthly_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = XGBOOST,
            formula = f,

            nrounds = 100,
            max_depth = 6,
            learning_rate = 0.001,
            min_child_weight = 1,
            subsample = 0.2,
            colsample_bytree = 1,

            idp = 1.5,
            nmax = 15
        )$var1.pred

    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    monthly_rain_plot(
        grid_geodata,
        monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/historic/IDW_XGB/IDW_XGB_",
            current_year, "_",
            sprintf("%02d", current_month),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 2] <- et - st

    print(date_index)
}

print("Mean Times")
print(apply(times, 2, mean))
