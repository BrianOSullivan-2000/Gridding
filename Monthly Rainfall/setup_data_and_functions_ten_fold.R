
## Load earlier dependencies
load("Data/Monthly_Rainfall/train_test_80_20_2016-2025.RData")


## Setup for the five-fold data

rain_data <- read.csv(
    "Data/Monthly_Rainfall/folds/10fold_2016-2025.csv"
) |>
    left_join(station_geodata, by = "stno") |>
    left_join(station_LTAs_9120, by = "stno") |>
    mutate(
        predicted_rain = NA
    )

dates <- rain_data |>
    group_by(year, month) |>
    dplyr::select(year, month) |>
    slice(1)


get_monthly_rain_data <- function(
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
        dplyr::select(year, month) |>
        dplyr::distinct()

    monthly_rain_data <- rain_data |>
        semi_join(selected_dates, by = c("year", "month")) |>
        rowwise() |>
        mutate(
            LTA = get(paste0("m_", month)),
            y = rain / LTA
        ) |>
        ungroup()

    if (experiment_type == "train_test") {
        if (is.null(current_fold)) {
            stop("Need to give fold for train-test split")
        }

        monthly_rain_data <- list(
            train = monthly_rain_data |>
                filter(fold != current_fold),
            test = monthly_rain_data |>
                filter(fold == current_fold)
        )
    }

    monthly_rain_data
}




update_monthly_predictions <- function(
    rain_data,
    y_hat,
    dates,
    date_index = NULL,
    current_fold,
    date_indices = date_index
) {
    monthly_rain_data <- get_monthly_rain_data(
        rain_data,
        dates,
        current_fold = current_fold,
        date_indices = date_indices
    )

    test_data <- monthly_rain_data$test |>
        mutate(
            predicted_rain = y_hat * LTA
        )

    selected_dates <- dates[date_indices, , drop = FALSE] |>
        dplyr::select(year, month) |>
        dplyr::distinct()

    target_rows <- rain_data |>
        mutate(.row_id = row_number()) |>
        filter(fold == current_fold) |>
        semi_join(selected_dates, by = c("year", "month")) |>
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

# save.image(file = "Data/Monthly_Rainfall/train_test_ten_fold_2016-2025.RData")