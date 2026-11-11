
## Validation tests for daily temperature grids       ##
## Random Forest spatial prediction framework (RFsp) ##

## Fit a regression model using Random Forests       ##
## Spatial effects are included in the model using   ##
## buffer distances to observations as predictors    ##

# %%

library(dplyr)
source("interpolation/RFsp.R")

## Load starting data
load("Data/Daily_Temperature/train_test_80_20_2025.RData")

# %%

## Regression Formula
## NOTE pass it as a string
base_formula <-
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            "exp25k",
            "n5", "s5", "e5", "w5",
            "aspect_10km",
            "slope_10km",
            sep = " + "
        )
    )
base_formula_maxt <- sub("^y", "y_maxt", base_formula)
base_formula_mint <- sub("^y", "y_mint", base_formula)

# %%

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        RFsp(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],
            base_formula = base_formula_maxt,
            mtry = 40,
            min.node.size = 2
        )$pred

    y_hat_mint <-
        RFsp(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],
            base_formula = base_formula_mint,
            mtry = 40,
            min.node.size = 2
        )$pred

    temperature_data <-
        update_daily_predictions(
            temperature_data,
            y_hat_maxt,
            y_hat_mint,
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
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
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

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data,
            dates,
            date_index,
            experiment_type = "all_data"
        )

    y_hat_maxt <-
        RFsp(
            daily_temperature_data,
            grid_geodata,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            base_formula = base_formula_maxt,
            mtry = 40,
            min.node.size = 2
        )$pred

    y_hat_mint <-
        RFsp(
            daily_temperature_data,
            grid_geodata,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            base_formula = base_formula_mint,
            mtry = 40,
            min.node.size = 2
        )$pred

    grid_geodata$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    grid_geodata$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        grid_geodata,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/RFsp/RFsp_",
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
