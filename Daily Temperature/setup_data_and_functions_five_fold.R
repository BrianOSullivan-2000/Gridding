
## Load earlier dependencies
load("Data/Daily_Temperature/train_test_80_20_2025.RData")


## Setup for the five-fold data

temperature_data <- read.csv(
    "Data/Daily_Temperature/folds/ST_5fold_2025.csv"
) |>
    left_join(station_geodata, by = "stno") |>
    left_join(station_LTAs_max_9120, by = "stno") |>
    dplyr::select(-ann) |>
    dplyr::rename_with(
        \(x) paste0("maxt_", x),
        starts_with("m_")
    ) |>
    left_join(station_LTAs_min_9120, by = "stno") |>
    dplyr::select(-ann) |>
    dplyr::rename_with(
        \(x) paste0("mint_", x),
        starts_with("m_")
    ) |>
    mutate(
        predicted_maxt = NA,
        predicted_mint = NA
    )

dates <- temperature_data |>
    group_by(year, month, day) |>
    dplyr::select(year, month, day) |>
    slice(1)




get_daily_temperature_data <- function(
    temperature_data,
    dates,
    date_index = NULL,
    current_fold = NULL,
    experiment_type = c("train_test", "all_data"),
    date_indices = date_index
) {
    experiment_type <- match.arg(experiment_type)

    if (is.null(date_indices) || !length(date_indices)) {
        stop("Provide one or more date_indices")
    }
    if (anyNA(date_indices) ||
            any(date_indices < 1 | date_indices > nrow(dates))) {
        stop("date_indices contains an invalid index")
    }

    selected_dates <- dates[date_indices, , drop = FALSE] |>
        dplyr::select(year, month, day) |>
        dplyr::distinct()

    daily_temperature_data <- temperature_data |>
        semi_join(selected_dates, by = c("year", "month", "day")) |>
        rowwise() |>
        mutate(
            LTA_maxt = get(paste0("maxt_m_", month)),
            LTA_mint = get(paste0("mint_m_", month)),
            y_maxt = maxt - LTA_maxt,
            y_mint = mint - LTA_mint
        ) |>
        ungroup()

    if (experiment_type == "train_test") {
        if (is.null(current_fold)) {
            stop("Need to give fold for train-test split")
        }

        daily_temperature_data <- list(
            train = daily_temperature_data |>
                filter(fold != current_fold),
            test = daily_temperature_data |>
                filter(fold == current_fold)
        )
    }

    daily_temperature_data
}




update_daily_predictions <- function(
    temperature_data,
    y_hat_maxt,
    y_hat_mint,
    dates,
    date_index = NULL,
    current_fold,
    date_indices = date_index
) {
    daily_temperature_data <- get_daily_temperature_data(
        temperature_data,
        dates,
        current_fold = current_fold,
        date_indices = date_indices
    )

    test_data <- daily_temperature_data$test |>
        mutate(
            predicted_maxt = y_hat_maxt + LTA_maxt,
            predicted_mint = y_hat_mint + LTA_mint
        )

    selected_dates <- dates[date_indices, , drop = FALSE] |>
        dplyr::select(year, month, day) |>
        dplyr::distinct()

    target_rows <- temperature_data |>
        mutate(.row_id = row_number()) |>
        filter(fold == current_fold) |>
        semi_join(selected_dates, by = c("year", "month", "day")) |>
        pull(.row_id)

    if (length(target_rows) != nrow(test_data)) {
        stop("Prediction count does not match the selected test rows")
    }

    temperature_data$predicted_maxt[target_rows] <- test_data$predicted_maxt
    temperature_data$predicted_mint[target_rows] <- test_data$predicted_mint

    temperature_data
}




get_temperature_data <- function(
    temperature_data,
    current_fold = NULL,
    experiment_type = c("train_test", "all_data"),
    date_index = NULL,
    dates = NULL,
    date_indices = date_index
) {
    experiment_type <- match.arg(experiment_type)

    if (!is.null(date_indices)) {
        if (is.null(dates)) {
            stop("Provide dates when using date_indices")
        }
        if (anyNA(date_indices) ||
                any(date_indices < 1 | date_indices > nrow(dates))) {
            stop("date_indices contains an invalid index")
        }

        selected_dates <- dates[date_indices, , drop = FALSE] |>
            dplyr::select(year, month, day) |>
            dplyr::distinct()
    }

    prepared_temperature_data <- temperature_data |>
        rowwise() |>
        mutate(
            LTA_maxt = get(paste0("maxt_m_", month)),
            LTA_mint = get(paste0("mint_m_", month)),
            y_maxt = maxt - LTA_maxt,
            y_mint = mint - LTA_mint
        ) |>
        ungroup() |>
        group_by(year, month, day) |>
        mutate(t = cur_group_id()) |>
        ungroup()

    if (experiment_type == "train_test") {
        if (is.null(current_fold)) {
            stop("Need to give fold for train-test split")
        }

        prepared_temperature_data <- list(
            train = prepared_temperature_data |>
                filter(fold != current_fold),
            test = prepared_temperature_data |>
                filter(fold == current_fold)
        )

        if (!is.null(date_indices)) {
            prepared_temperature_data$test <-
                prepared_temperature_data$test |>
                semi_join(selected_dates, by = c("year", "month", "day"))
        }
    } else if (!is.null(date_indices)) {
        prepared_temperature_data <- prepared_temperature_data |>
            semi_join(selected_dates, by = c("year", "month", "day"))
    }

    prepared_temperature_data
}




update_predictions <- function(
    temperature_data,
    y_hat_maxt,
    y_hat_mint,
    current_fold,
    dates = NULL,
    date_indices = NULL
) {
    prepared_temperature_data <- get_temperature_data(
        temperature_data,
        current_fold = current_fold,
        dates = dates,
        date_indices = date_indices
    )

    test_data <- prepared_temperature_data$test |>
        mutate(
            predicted_maxt = y_hat_maxt + LTA_maxt,
            predicted_mint = y_hat_mint + LTA_mint
        )

    target_data <- temperature_data |>
        mutate(.row_id = row_number()) |>
        filter(fold == current_fold)

    if (!is.null(date_indices)) {
        selected_dates <- dates[date_indices, , drop = FALSE] |>
            dplyr::select(year, month, day) |>
            dplyr::distinct()

        target_data <- target_data |>
            semi_join(selected_dates, by = c("year", "month", "day"))
    }

    target_rows <- target_data$.row_id

    if (length(target_rows) != nrow(test_data)) {
        stop("Prediction count does not match the selected test rows")
    }

    temperature_data$predicted_maxt[target_rows] <- test_data$predicted_maxt
    temperature_data$predicted_mint[target_rows] <- test_data$predicted_mint

    temperature_data
}




collect_metrics <- function(
    metrics_table,
    name,
    temperature_data
) {
    require(Metrics)

    rows <- list()
    for (element_idx in 1:2) {
        Element_name <-
            c("Min Temperature", "Max Temperature")[element_idx]
        element <- c("mint", "maxt")[element_idx]
        RMSEs <- c()
        MEs <- c()
        MAEs <- c()
        R2s <- c()
        ME10s <- c()

        for (current_fold in unique(temperature_data$fold)) {
            y <- temperature_data[
                temperature_data$fold == current_fold,
                element
            ]
            y_hat <- temperature_data[
                temperature_data$fold == current_fold,
                paste0("predicted_", element)
            ]

            RMSEs <- c(RMSEs, rmse(y_hat, y))
            MEs <- c(MEs, me(y_hat, y))
            MAEs <- c(MAEs, mae(y_hat, y))
            R2s <- c(R2s, R2(y_hat, y))
            ME10s <- c(
                ME10s,
                as.numeric(strsplit(ME10(y_hat, y), ",")[[1]])
            )
        }

        rows[[element_idx]] <- data.frame(
            Element = Element_name,
            Name = name,
            RMSE = mean(RMSEs),
            ME = mean(MEs),
            MAE = mean(MAEs),
            R2 = mean(R2s),
            ME10 = paste(
                sort(ME10s, decreasing = TRUE)[1:10],
                collapse = ", "
            )
        )
    }

    if (missing(metrics_table) || is.null(metrics_table)) {
        bind_rows(rows)
    } else {
        bind_rows(
            metrics_table,
            rows[[1]],
            rows[[2]]
        )
    }
}




## This will help make sure the earlier dependencies
## are up to date
cat(
    "setup_data_and_functions last updated on",
    last_update, "\n"
)

# save.image(file = "Data/Daily_Temperature/train_test_five_fold_2025.RData")