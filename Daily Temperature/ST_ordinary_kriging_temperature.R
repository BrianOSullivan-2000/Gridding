
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting ##
## Simple method that uses weighted average based on inverse distance ##

## For spatio-temporal data, consider time as a third dimension, scaled
## using the hhyperparameter C

# %%

library(dplyr)
source("interpolation/spatio_temporal_kriging.R")
source("interpolation/ordinary_kriging.R")
source("utilities/variogram.R")

## Load starting data
load("Data/Daily_Temperature/train_test_five_fold_2025.RData")

# %%

## First just spatial kriging
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
            ordinary_kriging(
                daily_temperature_data$train$y_maxt,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                vgm_model = "Sph",
                cutoff = 400000,
                width = 20000,
                nmax = Inf,
                debug.level = 0,
                plot_vgm = TRUE
            )$var1.pred

        y_hat_mint <-
            ordinary_kriging(
                daily_temperature_data$train$y_mint,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                vgm_model = "Sph",
                cutoff = 400000,
                width = 20000,
                nmax = Inf,
                debug.level = 0
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
        "Ordinary Kriging Sph",
        temperature_data |>
            as.data.frame()
    )

print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/OK_Sph.csv",
#     row.names = FALSE
# )
