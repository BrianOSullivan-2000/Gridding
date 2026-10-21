
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

vgm_model <- "Mat"

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

            init_pars = c(0.2, 0.1, 300, 1),
            lower = c(0.05, 0, 50, 0.2),
            upper = c(0.8, 0.4, 800, 2.5)

        )$pred

    y_hat_mint <-
        MLE_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            cov_function = vgm_model,
            formula = f_mint,

            init_pars = c(0.2, 0.1, 300, 1),
            lower = c(0.05, 0, 50, 0.2),
            upper = c(0.8, 0.4, 800, 2.5)

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

write.csv(
    metrics_table,
    paste0(
        "Results/Daily_Temperature/train_test_80_20_2025/",
        "RK_", vgm_model, "_MLE_safe.csv"
    ),
    row.names = FALSE
)
