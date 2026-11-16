
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
        temperature_data |> as.data.frame()
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/IDW.csv",
#     row.names = FALSE
# )

# %%

## Now spatio-temporal IDW

## Note - I had initially checked a lot more hyperparameter combinations
## (but lost the results...). To save a bit of time, I'm doing a
## smaller grid search over the hyperparameters that did well previously
## (Just leaving this note to avoid doing redundant future tests)

## Set up grid search for hyperparameters
idps <- c(1, 1.5, 2)
nmaxs <- c(10, 15, 20)
Cs <- c(30000, 50000, 75000, 100000)
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
            temperature_data |>
                as.data.frame()
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

# %%

## Plotting Time!!

idp <- 1
nmax <- 15
C <- 75000

plot_grid <- grid_geodata

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    ## First spatial IDW
    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data,
            dates,
            date_index,
            experiment_type = "all_data"
        )

    y_hat_maxt <-
        IDW(
            daily_temperature_data$y_maxt,
            daily_temperature_data[c("east", "north")],
            plot_grid[c("east", "north")],
            idp = idp,
            nmax = nmax
        )$var1.pred

    y_hat_mint <-
        IDW(
            daily_temperature_data$y_mint,
            daily_temperature_data[c("east", "north")],
            plot_grid[c("east", "north")],
            idp = idp,
            nmax = nmax
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
            "IDW/IDW_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )


    ## ST IDW

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

    y_hat_maxt <-
        spatio_temporal_IDW(
            prepared_temperature_data$y_maxt,
            prepared_temperature_data[c("east", "north", "t")],
            plot_grid[c("east", "north", "t")],
            idp = idp,
            nmax = nmax,
            C = C
        )$var1.pred

    y_hat_mint <-
        spatio_temporal_IDW(
            prepared_temperature_data$y_mint,
            prepared_temperature_data[c("east", "north", "t")],
            plot_grid[c("east", "north", "t")],
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
            "ST_IDW/ST_IDW_",
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
