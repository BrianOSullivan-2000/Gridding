
## Validation tests for daily temperature grids ##

## Regression Kriging (MLE) ##

## Modelling rainfall as a sum of a linear trend with geographic
## covariates + a smooth Gaussian process with a spatial dependence

## This is the same as regression kriging, except both components
## are fit as a joint model using Maximum Likelihood Expectation

# %%

library(dplyr)
source("interpolation/MLE_kriging.R")

## Load starting data
load("Data/Daily_Temperature/train_test_80_20_2025.RData")

# %%

## Regression Formulae

f_maxt <- as.formula(
    paste(
        "y_maxt ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            # "aspect_10km",
            # "n5", "e5", "s5", "w5",
            # "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

f_mint <- as.formula(
    paste(
        "y_mint ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            # "aspect_10km",
            # "n5", "e5", "s5", "w5",
            # "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

# %%

vgm_model <- "Sph"

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        MLE_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            cov_function = vgm_model,
            formula = f_maxt,

            init_pars = c(0.2, 0.1, 300),
            lower = c(0.05, 0, 50),
            upper = c(0.8, 0.4, 800)

        )$pred

    y_hat_mint <-
        MLE_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            cov_function = vgm_model,
            formula = f_mint,

            init_pars = c(0.2, 0.1, 300),
            lower = c(0.05, 0, 50),
            upper = c(0.8, 0.4, 800)

        )$pred

    temperature_data <-
        update_daily_predictions(
            temperature_data,
            y_hat_maxt,
            y_hat_mint,
            dates,
            date_index
        )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        paste0(
            "MLE Kriging ", vgm_model
        ),
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "RK_", vgm_model, "_MLE.csv"
#     ),
#     row.names = FALSE
# )


# %%

## Plot grids and check computation time

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

vgm_model <- "Sph"

for (date_index in seq_len(nrow(dates))[32:200]) {

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
        MLE_kriging(
            daily_temperature_data,
            grid_geodata,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            cov_function = vgm_model,
            formula = f_maxt,

            init_pars = c(0.2, 0.1, 300),
            lower = c(0.05, 0, 50),
            upper = c(0.8, 0.4, 800)
        )$pred

    y_hat_mint <-
        MLE_kriging(
            daily_temperature_data,
            grid_geodata,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            cov_function = vgm_model,
            formula = f_mint,

            init_pars = c(0.2, 0.1, 300),
            lower = c(0.05, 0, 50),
            upper = c(0.8, 0.4, 800)
        )$pred


    plot_grid$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    plot_grid$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        plot_grid,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/",
            "MLE_Kriging_", vgm_model,
            "/MLE_Kriging_", vgm_model, "_",
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
        mean(times)
    )
)
