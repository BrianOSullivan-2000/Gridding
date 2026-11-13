
## Validation tests for daily rainfall grids       ##
## Random Forest spatial prediction framework (RFsp) ##

## Fit a regression model using Random Forests       ##
## Spatial effects are included in the model using   ##
## buffer distances to observations as predictors    ##

# %%

library(dplyr)
source("interpolation/RFsp.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_provisional.RData")

# %%

## Regression Formula
## NOTE pass it as a string
base_formula <-
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            # "dist2c",
            # "exp25k",
            # "n5", "s5", "e5", "w5",
            # "aspect_10km_cos", "aspect_10km_sin",
            sep = " + "
        )
    )

# %%

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <-
        get_daily_rain_data(
            rain_data, dates, date_index
        )

    y_hat <-
        RFsp(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            base_formula = base_formula,
            mtry = 20,
            min.node.size = 2
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
        "RFsp",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_provisional/",
#         "RFsp.csv"
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
        RFsp(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            base_formula = base_formula,
            mtry = 20,
            min.node.size = 2
        )$pred

    grid_geodata$rain <-
        y_hat * grid_LTAs_9120[[paste0("m_", current_month)]]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall_provisional/RFsp/RFsp_",
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
