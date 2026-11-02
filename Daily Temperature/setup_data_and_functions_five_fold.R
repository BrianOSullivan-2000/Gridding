
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
    date_index,
    current_fold = NULL,
    experiment_type = c("train_test", "all_data")
) {

    experiment_type <- match.arg(experiment_type)

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    daily_temperature_data <- temperature_data |>
        filter(
            day == current_day,
            month == current_month,
            year == current_year
        ) |>
        mutate(
            LTA_maxt = .data[[paste0("maxt_m_", current_month)]],
            LTA_mint = .data[[paste0("mint_m_", current_month)]],
            y_maxt = maxt - LTA_maxt,
            y_mint = mint - LTA_mint
        )

    if (experiment_type == "train_test") {
        if (missing(current_fold)) {
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
    date_index,
    current_fold
) {
    daily_temperature_data <- get_daily_temperature_data(
        temperature_data,
        dates, date_index,
        current_fold = current_fold
    )

    daily_temperature_data$test <- daily_temperature_data$test |>
        mutate(
            predicted_maxt = y_hat_maxt + LTA_maxt,
            predicted_mint = y_hat_mint + LTA_mint
        )

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    temperature_data[
        temperature_data$day == current_day &
            temperature_data$month == current_month &
            temperature_data$year == current_year &
            temperature_data$fold == current_fold,

        "predicted_maxt"
    ] <-
        daily_temperature_data$test$predicted_maxt

    temperature_data[
        temperature_data$day == current_day &
            temperature_data$month == current_month &
            temperature_data$year == current_year &
            temperature_data$fold == current_fold,

        "predicted_mint"
    ] <-
        daily_temperature_data$test$predicted_mint

    temperature_data
}




get_temperature_data <- function(
    temperature_data,
    current_fold = NULL,
    experiment_type = c("train_test", "all_data"),
    date_index = NULL,
    dates = NULL
) {

    experiment_type <- match.arg(experiment_type)

    if (!missing(date_index) && !missing(dates)) {

        current_day <- dates[date_index, ]$day
        current_month <- dates[date_index, ]$month
        current_year <- dates[date_index, ]$year

        prepared_temperature_data <- temperature_data |>
            mutate(
                LTA_maxt = .data[[paste0("maxt_m_", current_month)]],
                LTA_mint = .data[[paste0("mint_m_", current_month)]]
            )

    } else {

        prepared_temperature_data <- temperature_data |>
            rowwise() |>
            mutate(
                LTA_maxt = get(paste0("maxt_m_", month)),
                LTA_mint = get(paste0("mint_m_", month))
            ) |>
            ungroup()

    }

    prepared_temperature_data <- prepared_temperature_data |>
        mutate(
            y_maxt = maxt - LTA_maxt,
            y_mint = mint - LTA_mint
        ) |>
        group_by(year, month, day) |>
        mutate(t = cur_group_id()) |>
        ungroup()

    if (experiment_type == "train_test") {
        if (missing(current_fold)) {
            stop("Need to give fold for train-test split")
        }

        prepared_temperature_data <- list(
            train = prepared_temperature_data |>
                filter(fold != current_fold),
            test = prepared_temperature_data |>
                filter(fold == current_fold)
        )
    }

    if (!missing(date_index) && !missing(dates)) {
        prepared_temperature_data$test <-
            prepared_temperature_data$test |>
            filter(
                year == current_year,
                month == current_month,
                day == current_day
            )
    }

    prepared_temperature_data
}




update_predictions <- function(
    temperature_data,
    y_hat_maxt,
    y_hat_mint,
    current_fold
) {
    prepared_temperature_data <- get_temperature_data(
        temperature_data,
        current_fold = current_fold
    )

    prepared_temperature_data$test <- prepared_temperature_data$test |>
        mutate(
            predicted_maxt = y_hat_maxt + LTA_maxt,
            predicted_mint = y_hat_mint + LTA_mint
        )

    temperature_data[
        temperature_data$fold == current_fold,
        "predicted_maxt"
    ] <-
        prepared_temperature_data$test$predicted_maxt

    temperature_data[
        temperature_data$fold == current_fold,
        "predicted_mint"
    ] <-
        prepared_temperature_data$test$predicted_mint

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

save.image(file = "Data/Daily_Temperature/train_test_five_fold_2025.RData")