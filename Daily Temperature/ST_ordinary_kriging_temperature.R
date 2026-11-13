
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting ##
## Simple method that uses weighted average based on inverse distance ##

## For spatio-temporal data, consider time as a third dimension, scaled
## using the hhyperparameter C

# %%

library(dplyr)
source("interpolation/spatio_temporal_kriging.R")
source("utilities/variogram.R")

## Load starting data
load("Data/Daily_Temperature/train_test_five_fold_2025.RData")

# %%

## First just spatial IDW

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

    }
}
