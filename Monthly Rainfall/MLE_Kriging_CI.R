
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

for (date_index in seq_len(nrow(dates))) {

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

    print(paste0(
        "Check for negative - ",
        sum(CI_lower < 0) / nrow(grid_geodata)
    ))
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

# %%

set.seed(222)

## I also want to check if the held out values from CV
## fall within the CI at a 95% coverage

coverages <- rep(NA, nrow(dates))

for (date_index in seq_len(nrow(dates))) {

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    result <-
        MLE_kriging(
            df = monthly_rain_data$train,
            new_df = monthly_rain_data$test,
            coords = monthly_rain_data$train[c("east", "north")],
            new_coords = monthly_rain_data$test[c("east", "north")],
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

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index,
        CI_lower = CI_lower, CI_upper = CI_upper
    )

    print(paste("Date", date_index))

    rain_test <- rain_data$test |>
        filter(
            year == current_year,
            month == current_month
        )

    print(paste0(
        "RMSE: ", rmse(rain_test$rain, rain_test$predicted_rain)
    ))
    print(paste0(
        "Coverage: ",
        sum(
            (rain_test$rain < rain_test$pred_upper) &
                (rain_test$rain > rain_test$pred_lower)
        ) * 100 / nrow(rain_test)
    ))

    coverage <- sum(
        (rain_test$rain < rain_test$pred_upper) &
            (rain_test$rain > rain_test$pred_lower)
    ) * 100 / nrow(rain_test)

    coverages[date_index] <- coverage
}

## parametric bootstrap?
x <- rain_data$test |>
    group_by(year, month) |>
    summarise(n = n())

mean(coverages, na.rm = TRUE)
