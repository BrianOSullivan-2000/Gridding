
## Code used to check the spatio-temporal variogram
## for daily temperature. This is an initial check to see if
## there is potential to do ST-interpolation

# %%

library(dplyr)

## Load starting data
load("Data/Daily_Temperature/train_test_80_20_2025.RData")

# %%

## First let's just plot a spatial variogram
source("visualisation/plot_variogram.R")

## Some sample data
date_index <- 1
daily_temperature_data <-
    get_daily_temperature_data(
        temperature_data,
        dates,
        date_index,
        experiment_type = "all_data"
    )

plot_spatial_variogram(
    daily_temperature_data,
    response = "y_mint",
    plot_curve = TRUE,
    psill = 1,
    nugget = 3,
    range = 3e5,
    kappa = 1.5
)

## It works!

# %%

## On to the spatio-temporal variogram

## Need a complete time index to plot this
## So load in all the ROI temperature data
df <- read.csv("Data/Daily_Temperature/daily_temperature_ROI_NI_2025.csv")

## Get a time index and merge with east + north
df <- df |>
    group_by(year, month, day) |>
    mutate(t = cur_group_id()) |>
    left_join(
        station_geodata |>
            dplyr::select(east, north, stno),
        by = "stno"
    ) |>
    left_join(
        station_LTAs_max_9120,
        by = "stno"
    ) |>
    filter(month == 1) |>
    mutate(y = maxt - m_1)

plot_spatio_temporal_variogram(
    df,
    response = "maxt",
    time_label = "Days",
    temporal_cutoff = 24
)
