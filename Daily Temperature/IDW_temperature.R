
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting ##
## Simple method that uses weighted average based on inverse distance ##

# %%

library(dplyr)
source("interpolation/IDW.R")

## Load starting data
load("Data/Daily_Temperature/train_test_80_20_2025.RData")

# %%

## Set up grid search for hyperparameters
idps <- c(0.5, 0.75, 1, 1.5, 2, 2.5, 3)
nmaxs <- c(3, 5, 8, 10, 12, 15, 20, 25)
hyperparameters <- expand.grid(idp = idps, nmax = nmaxs)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    idp <- hyperparameters[hyperparameter_index, ]$idp
    nmax <- hyperparameters[hyperparameter_index, ]$nmax

    for (date_index in seq_len(nrow(dates))) {

        daily_temperature_data <-
            get_daily_temperature_data(
                temperature_data, dates, date_index
            )

        y_hat_maxt <-
            IDW(
                daily_temperature_data$train$y_maxt,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                idp = idp,
                nmax = nmax
            )$var1.pred

        y_hat_mint <-
            IDW(
                daily_temperature_data$train$y_mint,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                idp = idp,
                nmax = nmax
            )$var1.pred

        temperature_data <-
            update_daily_predictions(
                temperature_data,
                y_hat_maxt,
                y_hat_mint,
                dates,
                date_index
            )
    }

    metrics_table <-
        collect_metrics(
            metrics_table,
            paste0("IDW idp-", idp, " nmax-", nmax),
            temperature_data$test$maxt,
            temperature_data$test$predicted_maxt,
            temperature_data$test$mint,
            temperature_data$test$predicted_mint
        )
}

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_80_20_2025/IDW_config.csv",
#     row.names = FALSE
# )

# %%

## IDW that selects hyperparameters using CV
IDW_grid_search <- function(
    y, coords, new_coords,
    idps = c(1, 1.5, 2, 2.5, 3),
    nmaxs = c(3, 5, 8, 10, 12, 15, 20),
    test_ratio = 0.2, seed = 222
) {

    set.seed(seed)

    test_idx <- sample(seq_along(y), size = floor(test_ratio * length(y)))
    y_train <- y[-test_idx]
    y_test <- y[test_idx]

    configs <- expand.grid(idp = idps, nmax = nmaxs)
    rmses <- rep(NA, nrow(configs))

    for (i in seq_len(nrow(configs))) {
        idp <- configs[i, ]$idp
        nmax <- configs[i, ]$nmax

        y_hat <- IDW(
            y_train, coords[-test_idx, ], coords[test_idx, ],
            idp = idp, nmax = nmax
        )$var1.pred

        rmses[i] <- rmse(y_hat, y_train)
    }

    idp <- configs[which.min(rmses), ]$idp
    nmax <- configs[which.min(rmses), ]$nmax

    y_hat <- IDW(
        y, coords, new_coords,
        idp = idp, nmax = nmax
    )

    y_hat
}

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        IDW_grid_search(
            daily_temperature_data$train$y_maxt,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],
            idps = c(1, 1.5, 2),
            nmaxs = c(10, 12, 15, 20)
        )$var1.pred

    y_hat_mint <-
        IDW_grid_search(
            daily_temperature_data$train$y_mint,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],
            idps = c(1, 1.5, 2),
            nmaxs = c(10, 12, 15, 20)
        )$var1.pred

    temperature_data <-
        update_daily_predictions(
            temperature_data,
            y_hat_maxt,
            y_hat_mint,
            dates,
            date_index
        )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW Grid Search",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_80_20_2025/IDW_grid_search.csv",
#     row.names = FALSE
# )


# %%

## Plot grids and check computation time

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    st <- Sys.time()

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

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
            grid_geodata[c("east", "north")],
            idp = 1,
            nmax = 15
        )$var1.pred

    y_hat_mint <-
        IDW(
            daily_temperature_data$y_mint,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            idp = 1,
            nmax = 15
        )$var1.pred

    grid_geodata$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    grid_geodata$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        grid_geodata,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/IDW/AA_IDW_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)
}

print(
    paste(
        "Mean Time",
        mean(times)
    )
)
