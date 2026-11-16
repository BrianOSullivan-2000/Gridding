
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
        temperature_data |>
            as.data.frame()
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/IDW_RF.csv",
#     row.names = FALSE
# )

# %%

## Now spatio-temporal regression IDW

## Regression Formula
f_RF <- as.formula(
    paste(
        "y ~",
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
f_MLR <- as.formula((
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            sep = " + "
        )
    )
))

f_RF_maxt <- update(f_RF, y_maxt ~ .)
f_RF_mint <- update(f_RF, y_mint ~ .)
f_MLR_maxt <- update(f_MLR, y_maxt ~ .)
f_MLR_mint <- update(f_MLR, y_mint ~ .)

# %%

## Try two regression methods (MLR and Random Forect)
regression_method <- Random_Forest

if (all(deparse(regression_method) == deparse(MLR))) {
    f_maxt <- f_MLR_maxt
    f_mint <- f_MLR_mint
}
if (all(deparse(regression_method) == deparse(Random_Forest))) {
    f_maxt <- f_RF_maxt
    f_mint <- f_RF_mint
}

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
        "ST IDW RF",
        temperature_data |>
            as.data.frame()
    )

print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/ST_IDW_RF.csv",
#     row.names = FALSE
# )

# %%

## Experimental check
## What if we tried just one month at a time
## rather than the entire year all at once

## Regression Formula
f_XGB <- as.formula(
    paste(
        "y ~",
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

f_XGB_maxt <- update(f_XGB, y_maxt ~ .)
f_XGB_mint <- update(f_XGB, y_mint ~ .)

f_maxt <- f_XGB_maxt
f_mint <- f_XGB_mint

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
        temperature_data |>
            as.data.frame()
    )

print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/ST_IDW_XGB.csv",
#     row.names = FALSE
# )


# %%

## Let's make some plots

plot_grid <- grid_geodata

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    ## First spatial IDW RF
    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data,
            dates,
            date_index,
            experiment_type = "all_data"
        )

    y_hat_maxt <-
        regression_IDW(
            daily_temperature_data,
            plot_grid,
            daily_temperature_data[c("east", "north")],
            plot_grid[c("east", "north")],
            formula = spatial_f_maxt,
            regression_method = Random_Forest,
            mtry = 5,
            min.node.size = 1,
            idp = 1,
            nmax = 15
        )$var1.pred

    y_hat_mint <-
        regression_IDW(
            daily_temperature_data,
            plot_grid,
            daily_temperature_data[c("east", "north")],
            plot_grid[c("east", "north")],
            formula = spatial_f_mint,
            regression_method = Random_Forest,
            mtry = 5,
            min.node.size = 1,
            idp = 1,
            nmax = 15
        )$var1.pred

    plot_grid$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    plot_grid$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        plot_grid,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/five_fold_tests/",
            "IDW_RF/IDW_RF_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )


    ## ST IDW (MLR)
    f_maxt <- f_MLR_maxt
    f_mint <- f_MLR_mint

    prepared_temperature_data <-
        get_temperature_data(
            temperature_data,
            current_fold = fold_idx,
            experiment_type = "all_data"
        )

    plot_grid$t <- which(
        (dates$day == current_day) &
            (dates$month == current_month) &
            (dates$year == current_year)
    )

    idp <- 1
    nmax <- 15
    C <- 50000

    y_hat_maxt <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data,
            plot_grid,
            prepared_temperature_data[c("east", "north", "t")],
            plot_grid[c("east", "north", "t")],

            regression_method = MLR,
            formula = f_maxt,

            idp = idp,
            nmax = nmax,
            C = C
        )$var1.pred

    y_hat_mint <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data,
            plot_grid,
            prepared_temperature_data[c("east", "north", "t")],
            plot_grid[c("east", "north", "t")],

            regression_method = MLR,
            formula = f_mint,

            idp = idp,
            nmax = nmax,
            C = C
        )$var1.pred


    grid_geodata$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    grid_geodata$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        grid_geodata,
        prepared_temperature_data |>
            filter(month == current_month),
        plot_destination = paste0(
            "Figures/Daily_Temperature/five_fold_tests/",
            "ST_IDW_MLR/ST_IDW_MLR_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    ## ST IDW (RF)
    f_maxt <- f_RF_maxt
    f_mint <- f_RF_mint

    C <- 50000
    mtry <- 6
    min.node.size <- 1

    y_hat_maxt <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data |>
                filter(month == current_month),
            plot_grid,
            prepared_temperature_data |>
                filter(month == current_month) |>
                dplyr::select(east, north, t),
            plot_grid[c("east", "north", "t")],

            regression_method = Random_Forest,
            formula = update(f_maxt, . ~ . - month),

            idp = idp,
            nmax = nmax,
            C = C,
            mtry = mtry,
            min.node.size = min.node.size
        )$var1.pred

    y_hat_mint <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data |>
                filter(month == current_month),
            plot_grid,
            prepared_temperature_data |>
                filter(month == current_month) |>
                dplyr::select(east, north, t),
            plot_grid[c("east", "north", "t")],

            regression_method = Random_Forest,
            formula = update(f_mint, . ~ . - month),

            idp = idp,
            nmax = nmax,
            C = C,
            mtry = mtry,
            min.node.size = min.node.size
        )$var1.pred

    grid_geodata$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    grid_geodata$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        grid_geodata,
        prepared_temperature_data |>
            filter(month == current_month),
        plot_destination = paste0(
            "Figures/Daily_Temperature/five_fold_tests/",
            "ST_IDW_RF/ST_IDW_RF_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )


    ## ST IDW (XGB once per month)
    f_maxt <- f_XGB_maxt
    f_mint <- f_XGB_mint

    nrounds <- 100
    max_depth <- 6
    learning_rate <- 0.001
    min_child_weight <- 1
    subsample <- 0.2
    colsample_bytree <- 1

    y_hat_maxt <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data |>
                filter(month == current_month),
            plot_grid,
            prepared_temperature_data |>
                filter(month == current_month) |>
                dplyr::select(east, north, t),
            plot_grid[c("east", "north", "t")],
            regression_method = XGBOOST,
            formula = update(f_maxt, . ~ . - month),
            idp = idp,
            nmax = nmax,
            C = C,
            nrounds = nrounds,
            max_depth = max_depth,
            learning_rate = learning_rate,
            min_child_weight = min_child_weight,
            subsample = subsample,
            colsample_bytree = colsample_bytree
        )$var1.pred

    y_hat_mint <-
        spatio_temporal_regression_IDW(
            prepared_temperature_data |>
                filter(month == current_month),
            plot_grid,
            prepared_temperature_data |>
                filter(month == current_month) |>
                dplyr::select(east, north, t),
            plot_grid[c("east", "north", "t")],
            regression_method = XGBOOST,
            formula = update(f_mint, . ~ . - month),
            idp = idp,
            nmax = nmax,
            C = C,
            nrounds = nrounds,
            max_depth = max_depth,
            learning_rate = learning_rate,
            min_child_weight = min_child_weight,
            subsample = subsample,
            colsample_bytree = colsample_bytree
        )$var1.pred

    daily_temperature_plot(
        grid_geodata,
        daily_temperature_data |>
            filter(month == current_month),
        plot_destination = paste0(
            "Figures/Daily_Temperature/five_fold_tests/",
            "ST_IDW_XGB/ST_IDW_XGB_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    print(
        paste0(
            date_index, "    ",
            current_day, "-",
            current_month, "-",
            current_year
        )
    )
}
