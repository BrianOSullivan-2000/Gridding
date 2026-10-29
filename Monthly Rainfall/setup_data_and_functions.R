
## Code used to set up starting environment for each gridding method ##

rm(list = ls())
library(dplyr)

## Station geodata and LTAs
station_geodata <- read.csv(
    "Data/General/station_geodata.csv"
)
station_LTAs_6190 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv"
)
station_LTAs_9120 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv"
)


## Grid geodata and LTAs
grid_geodata <- read.csv(
    "Data/General/grid_geodata.csv"
)
grid_LTAs_6190 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1961-1990.csv"
)
grid_LTAs_9120 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1991-2020.csv"
)


## Island Outline for plotting
island_outline <- read.csv(
    "Data/General/Irlcoast.csv"
)


## 80/20 Test (2016-2025) ##

## Combine obs, geodata, and LTAs
train_data <- read.csv("Data/Monthly_Rainfall/folds/train_split_2016-2025.csv")
test_data <- read.csv("Data/Monthly_Rainfall/folds/test_split_2016-2025.csv")

rain_data <- list(
    train = train_data,
    test = test_data
) |>
    lapply(\(data) {
        data |>
            left_join(station_geodata, by = "stno") |>
            left_join(station_LTAs_9120, by = "stno") |>
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
    baseline = FALSE,
    CI_lower = NULL,
    CI_upper = NULL
) {
    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    if (!baseline) {
        monthly_rain_data$test <- monthly_rain_data$test |>
            mutate(predicted_rain = (y_hat) * LTA)

        if (!missing(CI_lower) && !missing(CI_upper)) {
            monthly_rain_data$test <- monthly_rain_data$test |>
                mutate(
                    pred_lower = CI_lower * LTA,
                    pred_upper = CI_upper * LTA
                )
        }
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

    if (!missing(CI_lower) && !missing(CI_upper) && !baseline) {
        rain_data$test[
            rain_data$test$month == current_month &
                rain_data$test$year == current_year,

            "pred_lower"
        ] <-
            monthly_rain_data$test$pred_lower
        rain_data$test[
            rain_data$test$month == current_month &
                rain_data$test$year == current_year,

            "pred_upper"
        ] <-
            monthly_rain_data$test$pred_upper
    }

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




monthly_rain_CI_plot <- function(
    grid,
    monthly_rain_data,
    plot_destination,
    CI_level = 0.95
) {
    require(paletteer)
    require(cowplot)

    rain_breaks <- c(
        0, 20, 50, 100, 150, 200, 250, 300, 400, 500, 600, 800
    )

    rain_plot <-
        base_grid_plot(
            grid,
            island_outline,
            response = "rain",
            stations = monthly_rain_data,
            station_size = 1.5,
            breaks = rain_breaks,
            bias = 0.5,
            title = "Monthly Rain (mm)"
        )


    grid$diff <- grid$rain - grid$CI_lower
    z <- qnorm(1 - (1 - 0.95) / 2)
    grid$sd <- grid$diff / z

    sd_plot <-
        base_grid_plot(
            grid,
            island_outline,
            response = "sd",
            breaks = c(seq(0, 70, by = 5)),
            colour_palette = paletteer_c("viridis::viridis", 10),

            title = "Standard Deviation (mm)"
        )


    grid$perc_diff <-
        (grid$rain - grid$CI_lower) * 100 / grid$rain

    perc_diff_plot <-
        base_grid_plot(
            grid,
            island_outline,
            response = "perc_diff",
            colour_palette = c(
                "#4989ff",
                "#72c0ff",
                "#f7f4f4",
                "#df9f9b",
                "#f03e32",
                "#920a00"
            ),
            breaks = c(seq(0, 50, by = 5)),

            title = "% Difference (Rain - SD)"
        )

    full_plot <-
        plot_grid(
            rain_plot,
            sd_plot,
            perc_diff_plot,
            ncol = 3
        )

    ggsave(
        plot_destination,
        full_plot,
        width = 21,
        height = 7
    )
}

# save.image(file = "Data/Monthly_Rainfall/train_test_80_20_2016-2025.RData")
