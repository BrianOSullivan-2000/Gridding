
## Inspecting daily rainfall grids ##
## This script prepares plots of all published grids           ##
## for the 200 test days specified in the 2016-2025 experiment ##

# %%

library(dplyr)

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

count <- 0

for (current_year in 2016:2025) {
  
  daily_rainfall_grids <- read.csv(
    paste0(
      "Data/Daily_Rainfall/Published_Grids/",
      "IRL_DLY_RR_", current_year,
      "_grid.csv"
    )
  ) 
  
  annual_dates <- dates |>
    filter(year == current_year)
  
  for (date_index in seq_len(nrow(annual_dates))) {

    current_day <- annual_dates[date_index, ]$day
    current_month <- annual_dates[date_index, ]$month
    
    daily_rain_data <- get_daily_rain_data(
      rain_data, annual_dates, date_index, experiment_type = "all_data"
    )
    
    daily_column <- paste0(
      "X", current_year,
      sprintf("%02d", current_month),
      sprintf("%02d", current_day)
    )
    
    plot_grid <- daily_rainfall_grids[c("east", "north", daily_column)]
    plot_grid$rain <- plot_grid[daily_column]
    
    daily_rain_plot(
      plot_grid,
      daily_rain_data,
      plot_destination = paste0(
        "Figures/Daily_Rainfall/Published_Grids/Met_",
        current_year, "_",
        sprintf("%02d", current_month), "_",
        sprintf("%02d", current_day),
        ".jpg"
      )
    )
  }
}