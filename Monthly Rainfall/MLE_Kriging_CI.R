
## Regression Kriging (MLE) ##

## For Regression Kriging fit using Maximum Likelihood Estimation
## i.e. MLE_Kriging_Rainfall.R, you can get confidence intervals

## This notebook investigates plotting them and seeing if
## the held-out values fall within the 95% range during CV

# %%

library(dplyr)

source("interpolation/MLE_kriging.R")

## Load starting data
load("Data/Monthly_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## Regression Formula
f <- as.formula(
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            # "dist2c", "exp25k",
            sep = " + "
        )
    )
)

# %%

## Check grids and computation time


plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))[1:10]) {

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    result <-
        MLE_kriging(
            df = monthly_rain_data,
            new_df = grid_geodata,
            coords = monthly_rain_data[c("east", "north")],
            new_coords = grid_geodata[c("east", "north")],
            cov_function = "Exp",
            formula = f,
            init_pars = c(0.001, 0.0001, 150),
            lower = c(1e-8, 0, 10),
            upper = c(0.05, 0.01, 500),
            return_CI = TRUE,
            CI_type = "trend",
            CI_level = 0.95
        )
    y_hat <- result$pred
    CI_lower <- result$CI_lower
    CI_upper <- result$CI_upper

    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[[paste0("m_", current_month)]]
    grid_geodata$CI_lower <-
        (CI_lower) * grid_LTAs_9120[[paste0("m_", current_month)]]
    grid_geodata$CI_upper <-
        (CI_upper) * grid_LTAs_9120[[paste0("m_", current_month)]]

    monthly_rain_CI_plot(
        grid = grid_geodata,
        monthly_rain_data = monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/",
            "Kriging_MLE_Exp_Confidence_Intervals/",
            "Kriging_MLE_Exp_CI_",
            current_year, "_",
            sprintf(fmt = "%02d", current_month),
            ".jpg"
        )
    )

    print(date_index)
}
