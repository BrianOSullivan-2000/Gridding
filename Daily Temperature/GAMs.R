
## Validation tests for daily temperature grids ##
## Generalised Additive Models ##

## Fit a regression model using Generalized Additive models ##
## Models temperature as a combination of additive smooths  ##
## depending on covariates

# %%

library(dplyr)
source("regression/GAM.R")
source("interpolation/regression_grid.R")

## Load starting data
load("Data/Daily_Temperature/train_test_80_20_2025.RData")

# %%

## Regression Formula
f_maxt <- as.formula(
    paste(
        "y_maxt ~",
        paste(
            "te(east, north, k=5)",
            "s(points5, k=3)",
            "dist2c",
            "s(aspect_10km, bs='cc', k=4)",
            # "s(slope_10km, k=4)",
            sep = " + "
        )
    )
)

f_mint <- as.formula(
    paste(
        "y_mint ~",
        paste(
            "te(east, north, k=5)",
            "s(points5, k=3)",
            "dist2c",
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

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_grid(
            daily_temperature_data$train,
            daily_temperature_data$test,
            regression_method = GAM,
            formula = f_maxt
        )$pred


    y_hat_mint <-
        regression_grid(
            daily_temperature_data$train,
            daily_temperature_data$test,
            regression_method = GAM,
            formula = f_mint
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
        "GAM",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "GAM.csv"
#     ),
#     row.names = FALSE
# )

print(metrics_table)

# %%

## The results are quite good (and fast!)
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
        regression_grid(
            daily_temperature_data,
            grid_geodata,
            regression_method = GAM,
            formula = f_maxt
        )$pred

    y_hat_mint <-
        regression_grid(
            daily_temperature_data,
            grid_geodata,
            regression_method = GAM,
            formula = f_mint
        )$pred

    grid_geodata$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    grid_geodata$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        grid_geodata,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/GAM/GAM_",
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
