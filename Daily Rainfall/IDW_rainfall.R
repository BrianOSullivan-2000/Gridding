
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


dates <- rain_data$train |>
    group_by(year, month, day) |>
    dplyr::select(year, month, day) |>
    slice(1)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    idp <- hyperparameters[hyperparameter_index, ]$idp
    nmax <- hyperparameters[hyperparameter_index, ]$nmax

    for (date_index in seq_len(nrow(dates))) {
        current_day <- dates[date_index, ]$day
        current_month <- dates[date_index, ]$month
        current_year <- dates[date_index, ]$year

        daily_rain_data <- list(
            train = rain_data$train,
            test = rain_data$test
        ) |>
            lapply(\(data) {
                data |>
                    filter(
                        day == current_day,
                        month == current_month,
                        year == current_year
                    ) |>
                    mutate(
                        LTA = .data[[paste0("m_", current_month)]],
                        normalized_rain = rain / LTA,
                        y = log1p(normalized_rain)
                    )
            })

        y_hat <-
            IDW(
                daily_rain_data$train$y,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],
                idp = idp,
                nmax = nmax
            )$var1.pred

        daily_rain_data$test <- daily_rain_data$test |>
            mutate(predicted_rain = expm1(y_hat) * LTA)

        rain_data$test[
            rain_data$test$day == current_day &
                rain_data$test$month == current_month &
                rain_data$test$year == current_year,

            "predicted_rain"
        ] <-
            daily_rain_data$test$predicted_rain
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
#     "Results/Daily_Rainfall/IDW_config_log.csv",
#     row.names = FALSE
# )

# %%

## Small extension, run a grid search to identify best hyperparameters

IDW_grid_search <- function(
    y, coords, new_coords,
    idps, nmaxs,
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
    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    daily_rain_data <- list(
        train = rain_data$train,
        test = rain_data$test
    ) |>
        lapply(\(data) {
            data |>
                filter(
                    day == current_day,
                    month == current_month,
                    year == current_year
                ) |>
                mutate(
                    LTA = .data[[paste0("m_", current_month)]],
                    normalized_rain = rain / LTA,
                    y = log1p(normalized_rain)
                )
        })

    y_hat <-
        IDW_grid_search(
            daily_rain_data$train$y,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            idps = c(1.5, 2, 2.5),
            nmaxs = c(8, 10, 12, 15)
        )$var1.pred

    daily_rain_data$test <- daily_rain_data$test |>
        mutate(predicted_rain = expm1(y_hat) * LTA)

    rain_data$test[
        rain_data$test$day == current_day &
            rain_data$test$month == current_month &
            rain_data$test$year == current_year,

        "predicted_rain"
    ] <-
        daily_rain_data$test$predicted_rain
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
#     "Results/Daily_Rainfall/IDW_grid_search.csv",
#     row.names = FALSE
# )


# %%

## Plot grids and check computation time