
## Validation tests for daily rainfall grids ##

## Regression Kriging ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation with kriging, assuminig regression residuals   ##
##  follow a Gaussian Process (with spatial dependence)                   ##


# %%

library(dplyr)

source("interpolation/regression_kriging.R")

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
            "n5", "e5", "s5", "w5",
            "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

# %%

## Stepwise Regression

source("regression/MLR.R")

## Set up grid search for hyperparameters
cutoffs <- c(150000, 250000, 350000, 450000)
widths <- c(10000, 15000, 20000)
nmaxs <- c(8, 12, 15, 40, Inf)

hyperparameters <- expand.grid(cutoff = cutoffs, width = widths, nmax = nmaxs)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    cutoff <- hyperparameters[hyperparameter_index, ]$cutoff
    width <- hyperparameters[hyperparameter_index, ]$width
    nmax <- hyperparameters[hyperparameter_index, ]$nmax

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
                formula = f,

                cutoff = cutoff,
                width = width,
                nmax = nmax,

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
            paste0("RK cutoff-", cutoff, " width-", width, " nmax-", nmax),
            rain_data$test$rain,
            rain_data$test$predicted_rain
        )

    print(metrics_table)
}

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "RK_Exp_config_log.csv"
#     ),
#     row.names = FALSE
# )

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
            formula = f,

            cutoff = 350000,
            width = 15000,
            nmax = 15,

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
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "RK_Exp.csv"
#     ),
#     row.names = FALSE
# )

# %%

# Elastic-net Regularization

source("regression/elastic-net.R")

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <- tryCatch(
        regression_kriging(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = Elastic_Net,
            formula = f,

            cutoff = 350000,
            width = 15000,
            nmax = 15,

            alpha = 0.7,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred,

        error = function(e) {
            rep(0, length(daily_rain_data$test$y))
        }
    )

    rain_data <- update_daily_predictions(
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
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "RK_Exp_ENET.csv"
#     ),
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
        regression_kriging(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = GLS,
            formula = f,

            cutoff = 350000,
            width = 15000,
            nmax = 15,

            flex_fit = TRUE,
            vgm_model = "Sph",
            debug.level = 0
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "RK Sph GLS",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "RK_Sph_GLS_simple.csv"
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
            formula = f,

            nrounds = 300,
            max_depth = 6,
            learning_rate = 0.001,
            min_child_weight = 1,
            subsample = 0.2,
            colsample_bytree = 1,

            cutoff = 350000,
            width = 15000,
            nmax = 15,
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
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "RK_Exp_XGB_complex.csv"
#     ),
#     row.names = FALSE
# )
