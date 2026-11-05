
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting with regression ##

## For spatio-temporal data, consider time as a third dimension, scaled
## using the hhyperparameter C

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation, interpolate residuals using IDW ##

# %%

library(dplyr)
source("regression/MLR.R")
source("regression/random_forest.R")
source("regression/xgboost.R")
source("interpolation/regression_IDW.R")
source("interpolation/spatio_temporal_regression_IDW.R")

## Load starting data
load("Data/Daily_Temperature/train_test_five_fold_2025.RData")

# %%

## First spatial regression-IDW

## Regression Formula
spatial_f_maxt <- as.formula(
    paste(
        "y_maxt ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            "n5", "s5", "e5", "w5",
            "aspect_10km",
            sep = " + "
        )
    )
)

spatial_f_mint <- as.formula(
    paste(
        "y_mint ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            "n5", "s5", "e5", "w5",
            "aspect_10km",
            sep = " + "
        )
    )
)

for (fold_idx in seq_along(unique(temperature_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        # date_time <- Sys.time()

        daily_temperature_data <-
            get_daily_temperature_data(
                temperature_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat_maxt <-
            regression_IDW(
                daily_temperature_data$train,
                daily_temperature_data$test,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                formula = spatial_f_maxt,
                regression_method = Random_Forest,
                mtry = 5,
                min.node.size = 1,
                idp = 1,
                nmax = 15
            )$var1.pred

        y_hat_mint <-
            regression_IDW(
                daily_temperature_data$train,
                daily_temperature_data$test,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                formula = spatial_f_mint,
                regression_method = Random_Forest,
                mtry = 5,
                min.node.size = 1,
                idp = 1,
                nmax = 15
            )$var1.pred

        temperature_data <-
            update_daily_predictions(
                temperature_data,
                y_hat_maxt,
                y_hat_mint,
                dates,
                date_index,
                current_fold = fold_idx
            )

        # print(paste0(
        #     "Date ", date_index, " completed in ",
        #     Sys.time() - date_time
        # ))
    }

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}


metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW RF",
        temperature_data
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/IDW_RF.csv",
#     row.names = FALSE
# )

# %%

## Now spatio-temporal regression IDW

## Regression Formula
f_maxt <- as.formula(
    paste(
        "y_maxt ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            "n5", "s5", "e5", "w5",
            "aspect_10km",
            "month", "t",
            sep = " + "
        )
    )
)

f_mint <- as.formula(
    paste(
        "y_mint ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            "n5", "s5", "e5", "w5",
            "aspect_10km",
            "month", "t",
            sep = " + "
        )
    )
)

# %%

regression_method <- MLR
idp <- 1
nmax <- 15
C <- 75000
mtry <- 6
min.node.size <- 1

for (fold_idx in seq_along(unique(temperature_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    prepared_temperature_data <-
        get_temperature_data(
            temperature_data,
            current_fold = fold_idx
        )

    y_hat_maxt <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data$train,
            prepared_temperature_data$test,
            prepared_temperature_data$train[c("east", "north", "t")],
            prepared_temperature_data$test[c("east", "north", "t")],
            regression_method = regression_method,
            formula = f_maxt,
            idp = idp,
            nmax = nmax,
            C = C,
            mtry = mtry,
            min.node.size = min.node.size
        )$var1.pred

    y_hat_mint <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data$train,
            prepared_temperature_data$test,
            prepared_temperature_data$train[c("east", "north", "t")],
            prepared_temperature_data$test[c("east", "north", "t")],
            regression_method = regression_method,
            formula = f_maxt,
            idp = idp,
            nmax = nmax,
            C = C,
            mtry = mtry,
            min.node.size = min.node.size
        )$var1.pred

    temperature_data <-
        update_predictions(
            temperature_data,
            y_hat_maxt,
            y_hat_mint,
            current_fold = fold_idx
        )

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}


metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "ST IDW MLR",
        temperature_data
    )

print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/ST_IDW_MLR.csv",
#     row.names = FALSE
# )

# %%

## Experimental check
## What if we tried just one month at a time
## rather than the entire year all at once

regression_method <- XGBOOST
idp <- 1
nmax <- 15
C <- 50000

for (fold_idx in seq_along(unique(temperature_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (current_month in 1:12) {
        month_time <- Sys.time()

        date_indices <- which(dates$month == current_month)

        prepared_temperature_data <-
            get_temperature_data(
                temperature_data,
                current_fold = fold_idx,
                dates = dates,
                date_indices = date_indices
            )

        y_hat_maxt <-
            spatio_temporal_regression_IDW(
                prepared_temperature_data$train,
                prepared_temperature_data$test,
                prepared_temperature_data$train[c("east", "north", "t")],
                prepared_temperature_data$test[c("east", "north", "t")],
                regression_method = regression_method,
                formula = f_maxt,
                idp = idp,
                nmax = nmax,
                C = C,
                nrounds = 100,
                max_depth = 6,
                learning_rate = 0.001,
                min_child_weight = 1,
                subsample = 0.2,
                colsample_bytree = 1
            )$var1.pred

        y_hat_mint <-
            spatio_temporal_regression_IDW(
                prepared_temperature_data$train,
                prepared_temperature_data$test,
                prepared_temperature_data$train[c("east", "north", "t")],
                prepared_temperature_data$test[c("east", "north", "t")],
                regression_method = regression_method,
                formula = f_maxt,
                idp = idp,
                nmax = nmax,
                C = C,
                nrounds = 100,
                max_depth = 6,
                learning_rate = 0.001,
                min_child_weight = 1,
                subsample = 0.2,
                colsample_bytree = 1
            )$var1.pred

        temperature_data <-
            update_predictions(
                temperature_data,
                y_hat_maxt,
                y_hat_mint,
                current_fold = fold_idx,
                dates = dates,
                date_indices = date_indices
            )

        print(paste0(
            month.abb[current_month], " completed in ",
            Sys.time() - month_time
        ))
    }

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}


metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "ST IDW XGB",
        temperature_data
    )

print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/ST_IDW_XGB.csv",
#     row.names = FALSE
# )