
## Validation tests for daily rainfall grids ##
## Inverse Distance Weighting with regression ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation, interpolate residuals using IDW ##

# %%

library(dplyr)
source("interpolation/regression_IDW.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## Regression Formula
f <- as.formula(
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            "dist2c", "exp25k",
            "n5", "e5", "s5", "w5",
            "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

# %%

## Stepwise Regression

source("regression/MLR.R")

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = MLR,
            formula = f,
            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW MLR",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_MLR.csv",
#     row.names = FALSE
# )


# %%

# Elastic-net Regularization

source("regression/elastic-net.R")

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = Elastic_Net,
            formula = f,
            idp = 2,
            nmax = 15,
            alpha = 0.9
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW ENET",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_ENET.csv",
#     row.names = FALSE
# )


# %%

# Generalized Least Squares

source("regression/GLS.R")

for (date_index in seq_len(nrow(dates))) {
    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            regression_method = GLS,
            formula = f,
            flex_fit = TRUE, flex_vgm = TRUE,
            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW GLS",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_GLS.csv",
#     row.names = FALSE
# )

# %%

# Random Forests

source("interpolation/regression_IDW.R")
source("regression/random_forest.R")

for (date_index in seq_len(nrow(dates))) {
    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        regression_IDW(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],

            regression_method = Random_Forest,
            formula = f,

            mtry = c(5, 8, 10),
            min.node.size = c(2, 5, 10),
            cv_method = "oob",

            idp = 2,
            nmax = 15
        )$var1.pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
}

# write.csv(
#     metrics_table,
#     "Results/Daily_Rainfall/train_test_80_20_2016-2025/IDW_RF.csv",
#     row.names = FALSE
# )
