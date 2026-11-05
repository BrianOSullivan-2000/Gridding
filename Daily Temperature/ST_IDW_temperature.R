
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting ##
## Simple method that uses weighted average based on inverse distance ##

## For spatio-temporal data, consider time as a third dimension, scaled
## using the hhyperparameter C

# %%

library(dplyr)
source("interpolation/IDW.R")
source("interpolation/spatio_temporal_IDW.R")

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

        y_hat_maxt <-
            IDW(
                daily_temperature_data$train$y_maxt,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                idp = 1,
                nmax = 15
            )$var1.pred

        y_hat_mint <-
            IDW(
                daily_temperature_data$train$y_mint,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
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
        "IDW",
        temperature_data
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/IDW.csv",
#     row.names = FALSE
# )

# %%

## Now spatio-temporal IDW

## Set up grid search for hyperparameters
idps <- c(1, 1.5, 2)
nmaxs <- c(10, 15, 20, 30, 40, 100, 2000)
Cs <- c(5000, 10000, 20000, 30000, 40000, 50000, 75000, 100000)
hyperparameters <- expand.grid(idp = idps, nmax = nmaxs, C = Cs)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {
    print(paste0("Starting Hyperparameters ", hyperparameter_index))
    hyp_time <- Sys.time()

    idp <- hyperparameters[hyperparameter_index, ]$idp
    nmax <- hyperparameters[hyperparameter_index, ]$nmax
    C <- hyperparameters[hyperparameter_index, ]$C

    for (fold_idx in seq_along(unique(temperature_data$fold))) {

        print(paste0("Starting Fold ", fold_idx))
        fold_time <- Sys.time()

        prepared_temperature_data <-
            get_temperature_data(
                temperature_data,
                current_fold = fold_idx
            )

        y_hat_maxt <-
            spatio_temporal_IDW(
                prepared_temperature_data$train$y_maxt,
                prepared_temperature_data$train[c("east", "north", "t")],
                prepared_temperature_data$test[c("east", "north", "t")],
                idp = idp,
                nmax = nmax,
                C = C
            )$var1.pred

        y_hat_mint <-
            spatio_temporal_IDW(
                prepared_temperature_data$train$y_mint,
                prepared_temperature_data$train[c("east", "north", "t")],
                prepared_temperature_data$test[c("east", "north", "t")],
                idp = idp,
                nmax = nmax,
                C = C
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
            metrics_table,
            paste0("ST IDW idp-", idp, " nmax-", nmax, " C-", C),
            temperature_data
        )

    print(paste0(
        "Finished Hyperparameters ", hyperparameter_index,
        " in ", Sys.time() - hyp_time
    ))

}

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/ST_IDW_config.csv",
#     row.names = FALSE
# )
