
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting ##
## Simple method that uses weighted average based on inverse distance ##

# %%

library(dplyr)
source("interpolation/regression_kriging.R")

## Load starting data
load("Data/Daily_Temperature/train_test_80_20_2025.RData")

# %%

## Set up grid search for hyperparameters
## Set up grid search for hyperparameters
cutoffs <- c(200000, 250000, 300000, 350000, 400000)
widths <- c(5000, 10000, 15000, 20000, 25000)
hyperparameters <- expand.grid(cutoff = cutoffs, width = widths)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    cutoff <- hyperparameters[hyperparameter_index, ]$cutoff
    width <- hyperparameters[hyperparameter_index, ]$width

    for (date_index in seq_len(nrow(dates))) {

        daily_temperature_data <-
            get_daily_temperature_data(
                temperature_data, dates, date_index
            )

        vgm_model <- "Gau"

        y_hat_maxt <-
            ordinary_kriging(
                daily_temperature_data$train$y_maxt,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],

                cutoff = cutoff,
                width = width,
                flex_fit = TRUE,
                vgm_model = vgm_model,
                debug.level = 0
            )$var1.pred

        y_hat_mint <-
            ordinary_kriging(
                daily_temperature_data$train$y_mint,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],

                cutoff = cutoff,
                width = width,
                flex_fit = TRUE,
                vgm_model = vgm_model,
                debug.level = 0
            )$var1.pred

        temperature_data <-
            update_daily_predictions(
                temperature_data,
                y_hat_maxt,
                y_hat_mint,
                dates,
                date_index
            )
    }

    metrics_table <-
        collect_metrics(
            metrics_table,
            paste0("OK cutoff-", cutoff, " width-", width),
            temperature_data$test$maxt,
            temperature_data$test$predicted_maxt,
            temperature_data$test$mint,
            temperature_data$test$predicted_mint
        )

    print(hyperparameter_index)
}

write.csv(
    metrics_table,
    paste0(
        "Results/Daily_Temperature/train_test_80_20_2025/",
        "OK_Gau_config.csv"
    ),
    row.names = FALSE
)


# %%

## Plot grids and check computation time

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

    vgm_model <- "Sph"
    cutoff <- 400000
    width <- 20000
    nmax <- Inf

    y_hat_maxt <-
        ordinary_kriging(
            daily_temperature_data$y_maxt,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            cutoff = cutoff,
            width = width,
            flex_fit = TRUE,
            vgm_model = vgm_model,
            nmax = nmax,
            debug.level = 0
        )$var1.pred

    y_hat_maxt <-
        ordinary_kriging(
            daily_temperature_data$y_mint,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            cutoff = cutoff,
            width = width,
            flex_fit = TRUE,
            vgm_model = vgm_model,
            nmax = nmax,
            debug.level = 0
        )$var1.pred


    grid_geodata$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    grid_geodata$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        grid_geodata,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/",
            "Ordinary_Kriging_Sph/Ordinary_Kriging_Sph_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)
}

print(
    paste(
        "Mean Time",
        mean(times)
    )
)
