
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting with regression ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation, interpolate residuals using IDW ##

# %%

library(dplyr)
source("interpolation/regression_IDW.R")

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
            # "n5", "e5", "s5", "w5",
            # "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

# %%

## Stepwise Regression

source("regression/MLR.R")

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_IDW(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = MLR,
            formula = f_maxt,

            idp = idp,
            nmax = nmax
        )$var1.pred


    y_hat_mint <-
        regression_IDW(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = MLR,
            formula = f_mint,

            idp = 1,
            nmax = 15
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
        metrics_table = NULL,
        "IDW MLR",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "IDW_MLR.csv"
#     ),
#     row.names = FALSE
# )


# %%

# Elastic-net Regularization

source("regression/elastic-net.R")

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_IDW(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = Elastic_Net,
            formula = f_maxt,

            idp = idp,
            nmax = nmax,
            alpha = 0.1
        )$var1.pred


    y_hat_mint <-
        regression_IDW(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = Elastic_Net,
            formula = f_mint,

            idp = 1,
            nmax = 15,
            alpha = 0.1
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
        metrics_table = NULL,
        "IDW ENET",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "IDW_ENET.csv"
#     ),
#     row.names = FALSE
# )


# %%

# Generalised Least Squares

source("regression/GLS.R")

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_IDW(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = GLS,
            formula = f_maxt,
            flex_fit = TRUE, flex_vgm = TRUE,
            idp = 1,
            nmax = 15,

            vgm_model = "Sph"
        )$var1.pred


    y_hat_mint <-
        regression_IDW(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = GLS,
            formula = f_maxt,
            flex_fit = TRUE, flex_vgm = TRUE,
            idp = 1,
            nmax = 15,

            vgm_model = "Sph"
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
        metrics_table = NULL,
        "IDW GLS",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "IDW_GLS.csv"
#     ),
#     row.names = FALSE
# )
