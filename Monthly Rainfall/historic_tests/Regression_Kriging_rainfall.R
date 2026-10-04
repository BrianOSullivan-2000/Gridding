
## Validation tests for monthly rainfall grids ##

## Regression Kriging ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation with kriging, assuminig regression residuals   ##
##  follow a Gaussian Process (with spatial dependence)                   ##


# %%

library(dplyr)

source("interpolation/regression_kriging.R")

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
            # "dist2c", "exp25k",
            # "n5", "e5", "s5", "w5",
            # "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

# %%

source("regression/MLR.R")

for (date_index in seq_len(nrow(dates))) {

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_kriging(
            monthly_rain_data$train,
            monthly_rain_data$test,
            monthly_rain_data$train[c("east", "north")],
            monthly_rain_data$test[c("east", "north")],
            regression_method = MLR,
            formula = f,

            cutoff = 350000,
            width = 15000,
            nmax = 20,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred

    rain_data <- update_monthly_predictions(
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
#         "Results/Monthly_Rainfall/train_test_80_20_1883-1940/",
#         "RK_Exp.csv"
#     ),
#     row.names = FALSE
# )

# %%

# Elastic-net Regularization

source("regression/elastic-net.R")

for (date_index in seq_len(nrow(dates))) {

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    y_hat <- tryCatch(
        regression_kriging(
            monthly_rain_data$train,
            monthly_rain_data$test,
            monthly_rain_data$train[c("east", "north")],
            monthly_rain_data$test[c("east", "north")],
            regression_method = Elastic_Net,
            formula = f,

            cutoff = 350000,
            width = 15000,
            nmax = 20,

            alpha = 0.1,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred,

        error = function(e) {
            rep(0, length(monthly_rain_data$test$y))
        }
    )

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "RK Sph ENET",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Monthly_Rainfall/train_test_80_20_1883-1940/",
#         "RK_Exp_ENET.csv"
#     ),
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
        regression_kriging(
            monthly_rain_data$train,
            monthly_rain_data$test,
            monthly_rain_data$train[c("east", "north")],
            monthly_rain_data$test[c("east", "north")],

            regression_method = XGBOOST,
            formula = f,

            nrounds = 300,
            max_depth = 6,
            learning_rate = 0.001,
            min_child_weight = 1,
            subsample = 0.2,
            colsample_bytree = 1,

            cutoff = 350000,
            width = 15000,
            nmax = 20,
            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0

        )$var1.pred

    rain_data <- update_monthly_predictions(
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
#         "Results/Monthly_Rainfall/train_test_80_20_1883-1940/",
#         "RK_Exp_XGB.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Check grids and computation time

## Checking Exp with GLS (simple model)
## and Exp with XGBoost (simple model)

times <- matrix(nrow = 200, ncol = 2)
colnames(times) <- c("Exp MLR", "Exp XGB")

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    ## Exp MLR

    st <- Sys.time()

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_kriging(
            monthly_rain_data,
            grid_geodata,
            monthly_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = MLR,
            formula = y ~ east + north,

            cutoff = 350000,
            width = 15000,
            nmax = 20,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred

    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    monthly_rain_plot(
        grid_geodata,
        monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/historic/",
            "RK_MLR/RK_MLR_",
            current_year, "_",
            sprintf("%02d", current_month),
            ".jpg"
        )
    )

    et <- Sys.time()
    times[date_index, 1] <- et - st

    ## Exp XGBoost

    st <- Sys.time()

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        regression_kriging(
            monthly_rain_data,
            grid_geodata,
            monthly_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = XGBOOST,
            formula = y ~ east + north +
                points5 + dist2c + exp25k +
                n5 + e5 + s5 + w5 +
                ne5 + nw5 + se5 + sw5,

            nrounds = 300,
            max_depth = 6,
            learning_rate = 0.001,
            min_child_weight = 1,
            subsample = 0.2,
            colsample_bytree = 1,

            cutoff = 350000,
            width = 15000,
            nmax = 20,

            flex_fit = TRUE,
            vgm_model = "Exp",
            kappa = 0.5,
            debug.level = 0
        )$var1.pred

    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    monthly_rain_plot(
        grid_geodata,
        monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/historic/",
            "RK_Exp_XGB/RK_Exp_XGB_",
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
