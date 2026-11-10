
## Validation tests for daily rainfall grids ##
## Generalised Additive Models ##

## Fit a regression model using Generalized Additive models ##
## Models rainfall as a combination of additive smooths  ##
## depending on covariates

# %%

library(dplyr)
source("regression/GAM.R")
source("interpolation/regression_grid.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## Regression Formula
f <- as.formula(
    paste(
        "y ~",
        paste(
            "te(east, north, k=10)",
            "s(points5, k=5)",
            "dist2c",
            "ti(dist2c, east, k=5)",
            "s(aspect_10km, bs='cc', k=4)",
            # "s(slope_10km, k=4)",
            sep = " + "
        )
    )
)

# %%

## Generalized Additive Model

source("regression/GAM.R")

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <-
        get_daily_rain_data(
            rain_data, dates, date_index
        )

    y_hat <-
        regression_grid(
            daily_rain_data$train,
            daily_rain_data$test,
            regression_method = GAM,
            formula = f
        )$pred

    rain_data <-
        update_daily_predictions(
            rain_data,
            y_hat,
            dates,
            date_index
        )

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year
    print(
        paste0(
            date_index, "    ",
            current_day, "-",
            current_month, "-",
            current_year
        )
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "GAM",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "GAM.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Let's plot some grids and check computation time

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    st <- Sys.time()

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    daily_rain_data <-
        get_daily_rain_data(
            rain_data,
            dates,
            date_index,
            experiment_type = "all_data"
        )

    y_hat <-
        regression_grid(
            daily_rain_data,
            grid_geodata,
            regression_method = GAM,
            formula = f
        )$pred

    grid_geodata$rain <-
        y_hat * grid_LTAs_9120[[paste0("m_", current_month)]]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/GAM/GAM_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)

    print(
        paste0(
            date_index, "    ",
            current_day, "-",
            current_month, "-",
            current_year
        )
    )
}

print(
    paste(
        "Mean Time",
        mean(times, na.rm = TRUE)
    )
)
