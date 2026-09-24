
## Validation tests for daily rainfall grids ##
## Inverse Distance Weighting ##
## Simple method that uses weighted average based on inverse distance ##

# %%

library(dplyr)
source("interpolation/IDW.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## Set up grid search for hyperparameters
idps <- c(1, 1.5, 2, 2.5, 3)
nmaxs <- c(3, 5, 8, 10, 12, 15, 20)
hyperparameters <- expand.grid(idp = idps, nmax = nmaxs)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    idp <- hyperparameters[hyperparameter_index, ]$idp
    nmax <- hyperparameters[hyperparameter_index, ]$nmax

    for (date_index in seq_len(nrow(dates))) {

        daily_rain_data <- get_daily_rain_data(
            rain_data, dates, date_index
        )

        y_hat <-
            IDW(
                daily_rain_data$train$y,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],
                idp = idp,
                nmax = nmax
            )$var1.pred

        rain_data <- get_daily_predictions(
            rain_data, y_hat, dates, date_index
        )
    }

    metrics_table <-
        collect_metrics(
            metrics_table,
            paste0("IDW idp-", idp, " nmax-", nmax),
            rain_data$test$rain,
            rain_data$test$predicted_rain
        )
}

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_config_log.csv",
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

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        IDW_grid_search(
            daily_rain_data$train$y,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            idps = c(1.5, 2, 2.5),
            nmaxs = c(8, 10, 12, 15)
        )$var1.pred

    rain_data <- get_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW Grid Search",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_grid_search.csv",
#     row.names = FALSE
# )


# %%

## Plot grids and check computation time

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        IDW(
            daily_rain_data$y,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            idp = 2,
            nmax = 15
        )$var1.pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/IDW/IDW_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)
}

## I've added this to idp=2 nmax=15 row in the IDW results
print(paste(
    "Mean Time",
    mean(times)
))
