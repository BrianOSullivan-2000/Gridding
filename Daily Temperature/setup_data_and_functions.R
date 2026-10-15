
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
)
LTAs_max_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_max_grid_1991-2020.csv"
    )
)
LTAs_min_6190 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_grid_1961-1990.csv"
    )
)
LTAs_min_9120 <- read.csv(
    paste0(
        "Data/General/LTAs/Temperature/",
        "LTAs_Temperature_min_grid_1991-2020.csv"
    )
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
            new_row_maxt,
            new_row_mint
        )
    } else {
        bind_rows(
            metrics_table,
            new_row_maxt,
            new_row_mint
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


daily_temperature_plot <- function(
    grid, daily_temperature_data, plot_destination
) {

    require(paletteer)
    require(cowplot)

    full_breaks <- seq(-19, 34)
    full_palette <- colorRampPalette(
        colors = paletteer_c(
            palette = "grDevices::Zissou 1",
            n = length(full_breaks)
        )
    )(
        length(full_breaks)
    )

    max_breaks <- full_breaks[
        floor(min(grid$maxt)):ceiling(max(grid$maxt))
    ]
    min_breaks <- full_breaks[
        floor(min(grid$mint)):ceiling(max(grid$mint))
    ]

    max_palette <- full_palette[
        floor(min(grid$maxt)):ceiling(max(grid$maxt))
    ]
    min_palette <- full_palette[
        floor(min(grid$mint)):ceiling(max(grid$mint))
    ]

    maxt_plot <-
        base_grid_plot(
            grid,
            island_outline,
            response = "maxt",
            stations = daily_temperature_data,
            station_size = 1.5,
            breaks = max_breaks,
            colour_palette = max_palette
        )

    mint_plot <-
        base_grid_plot(
            grid,
            island_outline,
            response = "mint",
            stations = daily_temperature_data,
            station_size = 1.5,
            breaks = min_breaks,
            colour_palette = min_palette
        )

    full_plot <-
        plot_grid(maxt_plot, mint_plot, ncol = 2)

    ggsave(
        full_plot,
        plot_destination
    )
}

# Some palette options
# "grDevices::Zissou 1"
# "oompaBase::jetColors"
# "pals::warmcool"
# "pals::jet"
# "pals::kovesi.rainbow_bgyr_35_85_c73"



save.image(file = "Data/Daily_Temperature/train_test_80_20_2025.RData")
