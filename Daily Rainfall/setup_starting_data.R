
## Code used to set up starting environment for each gridding method ##

rm(list = ls())

station_geodata <- read.csv(
    "Data/General/station_geodata.csv"
)
station_LTAs_6190 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1961-1990.csv"
)
station_LTAs_9120 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_stations_1991-2020.csv"
)


grid_geodata <- read.csv(
    "Data/General/grid_geodata.csv"
)
grid_LTAs_6190 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1961-1990.csv"
)
grid_LTAs_9120 <- read.csv(
    "Data/General/LTAs/Rainfall/LTAs_Rainfall_grid_1991-2020.csv"
)

island_outline <- read.csv(
    "Data/General/Irlcoast.csv"
)


## 80/20 Test (2016-2025) ##

## Combine obs, geodata, and LTAs
train_data <- read.csv("Data/Daily_Rainfall/folds/train_split_2016-2025.csv")
test_data <- read.csv("Data/Daily_Rainfall/folds/test_split_2016-2025.csv")

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


source("visualisation/base_grid_plot.R")
source("visualisation/prediction_plot.R")
source("utilities/metric_functions.R")


metrics_table <- data.frame(
    Name = character(),
    RMSE = numeric(),
    ME = numeric(),
    Computation_Time = numeric(),
    R2 = numeric(),
    JSD = numeric(),
    ME10 = character(),
    Wet_Day_Bias = numeric()
)

collect_metrics <- function(metrics_table, name, y, y_hat) {
    new_row <- data.frame(
        Name = name,
        RMSE = rmse(y_hat, y),
        ME = me(y_hat, y),
        Computation_Time = NA_real_,
        MAE = mae(y_hat, y),
        R2 = R2(y_hat, y),
        JSD = JSD(y, y_hat),
        ME10 = ME10(y_hat, y),
        Wet_Day_Bias = Wet_Day_Bias(y_hat, y)
    )

    if (missing(metrics_table) || is.null(metrics_table)) {
        new_row
    } else {
        bind_rows(metrics_table, new_row)
    }
}

# save.image(file = "Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")
