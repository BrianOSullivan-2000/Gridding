
## Code used to set up starting environment for each gridding method ##

rm(list = ls())

## Station geodata and LTAs
station_geodata <- read.csv(
    "Data/General/station_geodata.csv"
)
station_LTAs_max_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_max_stations_1961-1990.csv"
    )
)
station_LTAs_max_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_max_stations_1991-2020.csv"
    )
)
station_LTAs_min_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_stations_1961-1990.csv"
    )
)
station_LTAs_min_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_stations_1991-2020.csv"
    )
)


## Grid geodata and LTAs
grid_geodata <- read.csv(
    "Data/General/grid_geodata.csv"
)
LTAs_max_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_max_grid_1961-1990.csv"
    )
) |>
    semi_join(
        grid_geodata,
        by = c("east", "north")
    )
LTAs_max_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_max_grid_1991-2020.csv"
    )
) |>
    semi_join(
        grid_geodata,
        by = c("east", "north")
    )
LTAs_min_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_grid_1961-1990.csv"
    )
) |>
    semi_join(
        grid_geodata,
        by = c("east", "north")
    )
LTAs_min_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_grid_1991-2020.csv"
    )
) |>
    semi_join(
        grid_geodata,
        by = c("east", "north")
    )


## Island Outline for plotting
island_outline <- read.csv(
    "Data/General/Irlcoast.csv"
)


## 80/20 Test (2016-2025) ##

## Combine obs, geodata, and LTAs
train_data <- read.csv("Data/Daily_Temperature/folds/train_split_2025.csv")
test_data <- read.csv("Data/Daily_Temperature/folds/test_split_2025.csv")

temperature_data <- list(
    train = train_data,
    test = test_data
) |>
    lapply(\(data) {
        data |>
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
    })

dates <- temperature_data$train |>
    group_by(year, month, day) |>
    dplyr::select(year, month, day) |>
    slice(1)


source("visualisation/base_grid_plot.R")
source("visualisation/prediction_plot.R")
source("utilities/metric_functions.R")



## Template for metrics
metrics_table <- data.frame(
    Element = character(),
    Name = character(),
    RMSE = numeric(),
    ME = numeric(),
    Computation_Time = numeric(),
    R2 = numeric(),
    ME10 = character()
)


## Function for getting metrics
collect_metrics <- function(
    metrics_table,
    name,
    y_maxt, y_hat_maxt,
    y_mint, y_hat_mint
) {

    require(Metrics)

    new_row_maxt <- data.frame(
        Element = "Max Temperature",
        Name = name,
        RMSE = rmse(y_hat_maxt, y_maxt),
        ME = me(y_hat_maxt, y_maxt),
        Computation_Time = NA_real_,
        MAE = mae(y_hat_maxt, y_maxt),
        R2 = R2(y_hat_maxt, y_maxt),
        ME10 = ME10(y_hat_maxt, y_maxt)
    )

    new_row_mint <- data.frame(
        Element = "Min Temperature",
        Name = name,
        RMSE = rmse(y_hat_mint, y_mint),
        ME = me(y_hat_mint, y_mint),
        Computation_Time = NA_real_,
        MAE = mae(y_hat_mint, y_mint),
        R2 = R2(y_hat_mint, y_mint),
        ME10 = ME10(y_hat_mint, y_mint)
    )

    if (missing(metrics_table) || is.null(metrics_table)) {
        rbind(
            new_row_mint,
            new_row_maxt
        )
    } else {
        bind_rows(
            metrics_table,
            new_row_mint,
            new_row_maxt
        )
    }
}


## Function that prepares data for a given day (for validation)
get_daily_temperature_data <- function(
    temperature_data,
    dates,
    date_index,
    experiment_type = c("train_test", "all_data")
) {

    experiment_type <- match.arg(experiment_type)

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    join_data <- function(data) {
        data <- data |>
            filter(
                day == current_day,
                month == current_month,
                year == current_year
            ) |>
            mutate(
                LTA_maxt = .data[[paste0("maxt_m_", current_month)]],
                LTA_mint = .data[[paste0("mint_m_", current_month)]],
                maxt_anomaly = maxt - LTA_maxt,
                mint_anomaly = mint - LTA_mint,
                y_maxt = maxt_anomaly,
                y_mint = mint_anomaly
            )
    }


    if (experiment_type == "train_test") {
        daily_temperature_data <- list(
            train = temperature_data$train,
            test = temperature_data$test
        )

        daily_temperature_data <- lapply(
            daily_temperature_data, join_data
        )

    } else if (experiment_type == "all_data") {
        daily_temperature_data <- bind_rows(
            temperature_data$train,
            temperature_data$test
        ) |>
            join_data()
    }

    daily_temperature_data
}


## Function that adds predicted temperature to test set
update_daily_predictions <- function(
    temperature_data,
    y_hat_maxt,
    y_hat_mint,
    dates,
    date_index
) {
    daily_temperature_data <- get_daily_temperature_data(
        temperature_data, dates, date_index
    )

    daily_temperature_data$test <- daily_temperature_data$test |>
        mutate(
            predicted_maxt = y_hat_maxt + LTA_maxt,
            predicted_mint = y_hat_mint + LTA_mint
        )

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    temperature_data$test[
        temperature_data$test$day == current_day &
            temperature_data$test$month == current_month &
            temperature_data$test$year == current_year,

        "predicted_maxt"
    ] <-
        daily_temperature_data$test$predicted_maxt

    temperature_data$test[
        temperature_data$test$day == current_day &
            temperature_data$test$month == current_month &
            temperature_data$test$year == current_year,

        "predicted_mint"
    ] <-
        daily_temperature_data$test$predicted_mint

    temperature_data
}


## Temperature palette pulled from swalsh scripts on sol
temperature_palette <- c(
    "#502020", "#611838", "#721050", "#830868", "#940080",
    "#000080", "#00009F", "#0000BF", "#0000DF", "#0000FF",
    "#566CF2", "#6781F0", "#7997ED", "#8AACEB", "#9BC2E8",
    "#91B7CE", "#87ACB5", "#7DA19B", "#739682", "#698B69",
    "#76AC76", "#7CBC7C", "#83CD83", "#89DD89", "#90EE90",
    "#A8D56A", "#B4C958", "#C1BD45", "#CDB132", "#DAA520",
    "#B1700D", "#C5740A", "#D87806", "#EB7C03", "#FF8000",
    "#FF5000", "#DF3C00", "#BF2800", "#9F1400", "#800000",
    "#993300", "#B26600", "#CC9900", "#E5CC00", "#FFFF00"
)
temperature_idx_sequence <- cbind(
    seq(-10, 34, 1), seq(1, 45, 1)
)


daily_temperature_plot <- function(
    grid, daily_temperature_data, plot_destination
) {

    require(cowplot)

    cut_maxt <- seq(
        trunc(min(grid$maxt)) - 1,
        trunc(max(grid$maxt)) + 1,
        1
    )
    cut_mint <- seq(
        trunc(min(grid$mint)) - 1,
        trunc(max(grid$mint)) + 1,
        1
    )

    sequence_maxt <- subset(
        temperature_idx_sequence,
        temperature_idx_sequence[, 1] >=
            min(cut_maxt) & temperature_idx_sequence[, 1] <= max(cut_maxt)
    )
    sequence_mint <- subset(
        temperature_idx_sequence,
        temperature_idx_sequence[, 1] >=
            min(cut_mint) & temperature_idx_sequence[, 1] <= max(cut_mint)
    )

    max_palette <- temperature_palette[
        min(sequence_maxt[, 2]):max(sequence_maxt[, 2])
    ]
    min_palette <- temperature_palette[
        min(sequence_mint[, 2]):max(sequence_mint[, 2])
    ]

    maxt_plot <-
        base_grid_plot(
            grid,
            island_outline,
            response = "maxt",
            stations = daily_temperature_data,
            grid_size = 0.07,
            station_size = 1.5,
            breaks = sequence_maxt[, 1],
            colour_palette = max_palette
        )

    mint_plot <-
        base_grid_plot(
            grid,
            island_outline,
            response = "mint",
            stations = daily_temperature_data,
            grid_size = 0.07,
            station_size = 1.5,
            breaks = sequence_mint[, 1],
            colour_palette = min_palette
        )

    full_plot <-
        plot_grid(maxt_plot, mint_plot, ncol = 2)

    ggsave(
        plot_destination,
        full_plot,
        width = 14,
        height = 7
    )

}

# save.image(file = "Data/Daily_Temperature/train_test_80_20_2025.RData")
