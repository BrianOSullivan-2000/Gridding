
## Code used to set up starting environment for each gridding method ##
## This time setting up for the historic monthly grids experiment    ##

rm(list = ls())
library(dplyr)

## Station geodata and LTAs
station_geodata <- bind_rows(
    read.csv(
        paste0(
            "Data/Monthly_Rainfall/historic_data/",
            "roi_stations_geo_data.csv"
        )
    ),
    read.csv(
        paste0(
            "Data/Monthly_Rainfall/historic_data/",
            "ni_stations_geo_data.csv"
        )
    )
) |>
    dplyr::rename(stno = STNO)

station_LTAs <- bind_rows(
    read.csv(
        paste0(
            "Data/Monthly_Rainfall/historic_data/",
            "roi_stations_LTA8110_data.csv"
        )
    ),
    read.csv(
        paste0(
            "Data/Monthly_Rainfall/historic_data/",
            "ni_stations_LTA8110_data.csv"
        )
    )
) |>
    dplyr::rename(stno = STNO) |>
    rename_with(
        ~ sub("rr8110m", "m_", .x),
        matches("^rr8110m")
    ) |>
    dplyr::select(-east, -north)


## Grid geodata and LTAs
grid_geodata <- read.csv(
    paste0(
        "Data/Monthly_Rainfall/historic_data/",
        "island_grid_geo_8110_ok.csv"
    )
) |>
    dplyr::select(
        east, north,
        niflag, elev, points5,
        n5, s5, e5, w5,
        ne5, se5, sw5, nw5,
        dist2c, exp5k, exp10k,
        exp15k, exp20k, exp25k
    )

grid_LTAs <- read.csv(
    paste0(
        "Data/Monthly_Rainfall/historic_data/",
        "island_grid_geo_8110_ok.csv"
    )
) |>
    rename_with(
        ~ sub("aar8110m", "m_", .x),
        matches("^aar8110m")
    ) |>
    dplyr::select(east, north, starts_with("m_"))


## Island Outline for plotting
island_outline <- read.csv(
    "Data/General/Irlcoast.csv"
)


## 80/20 Test (2016-2025) ##

## Combine obs, geodata, and LTAs
train_data <- read.csv("Data/Monthly_Rainfall/folds/train_split_1883-1940.csv")
test_data <- read.csv("Data/Monthly_Rainfall/folds/test_split_1883-1940.csv")

rain_data <- list(
    train = train_data,
    test = test_data
) |>
    lapply(\(data) {
        data |>
            left_join(station_geodata, by = "stno") |>
            left_join(station_LTAs, by = "stno") |>
            filter(if_all(everything(), ~ !is.na(.x))) |>
            mutate(predicted_rain = NA)
    })

dates <- rain_data$train |>
    group_by(year, month) |>
    dplyr::select(year, month) |>
    slice(1)


source("visualisation/base_grid_plot.R")
source("visualisation/prediction_plot.R")
source("utilities/metric_functions.R")


## Template for metrics
metrics_table <- data.frame(
    Name = character(),
    RMSE = numeric(),
    ME = numeric(),
    Computation_Time = numeric(),
    R2 = numeric(),
    JSD = numeric(),
    ME10 = character()
)


## Function for getting metrics
collect_metrics <- function(metrics_table, name, y, y_hat) {

    require(Metrics)

    new_row <- data.frame(
        Name = name,
        RMSE = rmse(y_hat, y),
        ME = me(y_hat, y),
        Computation_Time = NA_real_,
        MAE = mae(y_hat, y),
        R2 = R2(y_hat, y),
        JSD = JSD(y, y_hat),
        ME10 = ME10(y_hat, y)
    )

    if (missing(metrics_table) || is.null(metrics_table)) {
        new_row
    } else {
        bind_rows(metrics_table, new_row)
    }
}


## Function that prepares data for a given day (for validation)
get_monthly_rain_data <- function(
    rain_data,
    dates,
    date_index,
    experiment_type = c("train_test", "all_data")
) {

    experiment_type <- match.arg(experiment_type)

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    join_data <- function(data) {
        data <- data |>
            filter(
                month == current_month,
                year == current_year
            ) |>
            mutate(
                LTA = .data[[paste0("m_", current_month)]],
                normalized_rain = rain / LTA,
                y = (normalized_rain)
            )
    }


    if (experiment_type == "train_test") {
        monthly_rain_data <- list(
            train = rain_data$train,
            test = rain_data$test
        )

        monthly_rain_data <- lapply(
            monthly_rain_data, join_data
        )

    } else if (experiment_type == "all_data") {
        monthly_rain_data <- bind_rows(
            rain_data$train,
            rain_data$test
        ) |>
            join_data()
    }

    monthly_rain_data
}


## Function that takes y_hat (predicted rain) and adds it to test set
update_monthly_predictions <- function(
    rain_data,
    y_hat,
    dates,
    date_index,
    baseline = FALSE
) {
    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    if (!baseline) {
        monthly_rain_data$test <- monthly_rain_data$test |>
            mutate(predicted_rain = (y_hat) * LTA)
    } else {
        monthly_rain_data$test <- monthly_rain_data$test |>
            mutate(predicted_rain = y_hat)
    }

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    rain_data$test[
        rain_data$test$month == current_month &
            rain_data$test$year == current_year,

        "predicted_rain"
    ] <-
        monthly_rain_data$test$predicted_rain

    rain_data
}


monthly_rain_plot <- function(grid, monthly_rain_data, plot_destination) {

    base_grid_plot(
        grid,
        island_outline,
        response = "rain",

        stations = monthly_rain_data,
        station_size = 1.5,

        breaks = c(
            0, 20, 50, 100, 150, 200, 250, 300, 400, 500, 600, 800
        ),
        bias = 0.5,

        plot_destination = plot_destination
    )
}

# save.image(file = "Data/Monthly_Rainfall/train_test_80_20_1883-1940.RData")
