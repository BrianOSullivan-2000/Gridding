
## Climate normals were dervied from each day of the year from
## Met Éireann's monthly long-term averages for 1991-2020

## This script just plots them, if you want info on how they were
## generated, check out prepare_test_data.R

# %%

library(dplyr)

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

all_dates <- as.Date("2025-01-01") + 0:364

for (day_of_year in 1:365) {

    date_label <- format(all_dates[day_of_year], "%d %B")
    date_label <- sub("^0", "", date_label)

    daily_grid_path <- paste0(
        "Data/General/LTAs/Rainfall/Daily_Averages/",
        "grid_rainfall_9120_day", day_of_year, ".csv"
    )

    plot_grid <- grid_geodata[c("east", "north")]
    plot_grid$rain <- read.csv(daily_grid_path)[, 3]

    plot_destination <- paste0(
        "Figures/Daily_Rainfall/Normals/",
        "grid_rainfall_9120_day", day_of_year, ".jpg"
    )

    base_grid_plot(
        plot_grid,
        island_outline,
        response = "rain",

        station_size = 1.5,

        breaks = c(
            0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75,
            1, 2, 3, 4, 5, 7.5,
            10, 20, 30, 40, 50, 75, 100
        ),
        bias = 0.5,

        plot_destination = plot_destination,
        title = date_label,
        legend_title = "mm"
    )

    print(day_of_year)
}
