
## Load earlier dependencies
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")


## Setup for the five-fold data

rain_data <- read.csv(
    "Data/Daily_Rainfall/folds/10fold_2016-2025.csv"
) |>
    left_join(station_geodata, by = "stno") |>
    left_join(station_LTAs_9120, by = "stno") |>
    mutate(
        predicted_rain = NA
    )

dates <- rain_data |>
    group_by(year, month, day) |>
    dplyr::select(year, month, day) |>
    slice(1)


get_daily_rain_data <- function(
    rain_data,
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

    selected_dates <- dates[date_indices, , drop = FALSE] |>
        dplyr::select(year, month, day) |>
        dplyr::distinct()

    daily_rain_data <- rain_data |>
        semi_join(selected_dates, by = c("year", "month", "day")) |>
        rowwise() |>
        mutate(
            LTA = get(paste0("m_", month)),
            normalized_rain = rain / LTA,
            y = log1p(normalized_rain)
        ) |>
        ungroup()

    if (experiment_type == "train_test") {
        if (is.null(current_fold)) {
            stop("Need to give fold for train-test split")
        }

        daily_rain_data <- list(
            train = daily_rain_data |>
                filter(fold != current_fold),
            test = daily_rain_data |>
                filter(fold == current_fold)
        )
    }

    daily_rain_data
}




update_daily_predictions <- function(
    rain_data,
    y_hat,
    dates,
    date_index = NULL,
    current_fold,
    date_indices = date_index
) {
    daily_rain_data <- get_daily_rain_data(
        rain_data,
        dates,
        current_fold = current_fold,
        date_indices = date_indices
    )

    test_data <- daily_rain_data$test |>
        mutate(
            predicted_rain = expm1(y_hat) * LTA
        )

    selected_dates <- dates[date_indices, , drop = FALSE] |>
        dplyr::select(year, month, day) |>
        dplyr::distinct()

    target_rows <- rain_data |>
        mutate(.row_id = row_number()) |>
        filter(fold == current_fold) |>
        semi_join(selected_dates, by = c("year", "month", "day")) |>
        pull(.row_id)

    rain_data$predicted_rain[target_rows] <- test_data$predicted_rain

    rain_data
}




collect_metrics <- function(
    metrics_table,
    name,
    rain_data
) {
    require(Metrics)

    RMSEs <- c()
    MEs <- c()
    MAEs <- c()
    R2s <- c()
    Wet_Day_Biases <- c()
    ME10s <- c()

    for (current_fold in unique(rain_data$fold)) {
        y <- rain_data[
            rain_data$fold == current_fold,
            "rain"
        ]
        y_hat <- rain_data[
            rain_data$fold == current_fold,
            "predicted_rain"
        ]

        RMSEs <- c(RMSEs, rmse(y_hat, y))
        MEs <- c(MEs, me(y_hat, y))
        MAEs <- c(MAEs, mae(y_hat, y))
        R2s <- c(R2s, R2(y_hat, y))
        Wet_Day_Biases <- c(Wet_Day_Biases, Wet_Day_Bias(y_hat, y))
        ME10s <- c(
            ME10s,
            as.numeric(strsplit(ME10(y_hat, y), ",")[[1]])
        )
    }

    row <- data.frame(
        Name = name,
        RMSE = mean(RMSEs),
        ME = mean(MEs),
        MAE = mean(MAEs),
        R2 = mean(R2s),
        Wet_Day_Bias = mean(Wet_Day_Biases),
        ME10 = paste(
            sort(ME10s, decreasing = TRUE)[1:10],
            collapse = ", "
        )
    )

    if (missing(metrics_table) || is.null(metrics_table)) {
        row
    } else {
        bind_rows(
            metrics_table,
            row
        )
    }
}




## This will help make sure the earlier dependencies
## are up to date
cat(
    "setup_data_and_functions last updated on",
    last_update, "\n"
)

# save.image(file = "Data/Daily_Rainfall/train_test_ten_fold_2016-2025.RData")