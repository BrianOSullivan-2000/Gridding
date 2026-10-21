
## Validation tests for daily temperature grids ##

## Regression Kriging ##

## Step 1 is regression, fit trend to data based on geographic covariates ##
## Step 2 is interpolation with kriging, assuminig regression residuals   ##
##  follow a Gaussian Process (with spatial dependence)                   ##

# %%

library(dplyr)
source("interpolation/regression_kriging.R")

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
            "aspect_10km",
            "n5", "e5", "s5", "w5",
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
            "aspect_10km",
            "n5", "e5", "s5", "w5",
            # "ne5", "nw5", "se5", "sw5",
            sep = " + "
        )
    )
)

# %%

source("regression/MLR.R")

## Set up grid search for hyperparameters
cutoffs <- c(200000, 300000, 400000)
widths <- c(5000, 10000, 15000, 20000, 25000)
vgm_models <- c("Sph", "Exp", "Mat", "Gau")
hyperparameters <- expand.grid(
    cutoff = cutoffs, width = widths, vgm_model = vgm_models
)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    cutoff <- hyperparameters[hyperparameter_index, ]$cutoff
    width <- hyperparameters[hyperparameter_index, ]$width
    vgm_model <- as.character(
        hyperparameters[hyperparameter_index, ]$vgm_model
    )

    for (date_index in seq_len(nrow(dates))) {

        daily_temperature_data <-
            get_daily_temperature_data(
                temperature_data, dates, date_index
            )

        y_hat_maxt <-
            regression_kriging(
                daily_temperature_data$train,
                daily_temperature_data$test,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],

                regression_method = MLR,
                formula = f_maxt,

                cutoff = cutoff,
                width = width,
                vgm_model = vgm_model,

                nmax = Inf,
                flex_fit = TRUE,
                debug.level = 0
            )$var1.pred

        y_hat_mint <-
            regression_kriging(
                daily_temperature_data$train,
                daily_temperature_data$test,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],

                regression_method = MLR,
                formula = f_mint,

                cutoff = cutoff,
                width = width,
                vgm_model = vgm_model,

                nmax = Inf,
                flex_fit = TRUE,
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
            paste0(
                "OK cutoff-", cutoff,
                " width-", width,
                " vgm_model-", vgm_model
            ),
            temperature_data$test$maxt,
            temperature_data$test$predicted_maxt,
            temperature_data$test$mint,
            temperature_data$test$predicted_mint
        )

    print(hyperparameter_index)
}

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "RK_config.csv"
#     ),
#     row.names = FALSE
# )


# %%

source("regression/MLR.R")


cutoff <- 400000
width <- 20000
vgm_model <- "Exp"

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = MLR,
            formula = f_maxt,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
            debug.level = 0
        )$var1.pred

    y_hat_mint <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = MLR,
            formula = f_mint,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
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
        metrics_table = NULL,
        "Regression Kriging Exp",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "RK_Exp.csv"
#     ),
#     row.names = FALSE
# )

# %%

# Elastic-net Regularization

source("regression/elastic-net.R")

cutoff <- 400000
width <- 20000
alpha <- 0.1
vgm_model <- "Exp"

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = Elastic_Net,
            formula = f_maxt,

            alpha = alpha,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
            debug.level = 0
        )$var1.pred

    y_hat_mint <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = Elastic_Net,
            formula = f_mint,

            alpha = alpha,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
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
        metrics_table = NULL,
        "RK Exp ENET",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "RK_ENET.csv"
#     ),
#     row.names = FALSE
# )

# %%

# Generalised Least Squares

source("regression/GLS.R")

cutoff <- 400000
width <- 20000
vgm_model <- "Exp"

for (date_index in seq_len(nrow(dates))) {

    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = GLS,
            formula = f_maxt,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
            debug.level = 0
        )$var1.pred

    y_hat_mint <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = GLS,
            formula = f_mint,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
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
        metrics_table = NULL,
        "RK GLS",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "RK_GLS.csv"
#     ),
#     row.names = FALSE
# )

# %%

# Random Forests Grid Search

source("regression/random_forest.R")

cutoff <- 400000
width <- 20000
vgm_model <- "Exp"


for (date_index in seq_len(nrow(dates))) {
    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data, dates, date_index
        )

    y_hat_maxt <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = Random_Forest,
            formula = f_maxt,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
            debug.level = 0,

            mtry = 5,
            min.node.size = 1

        )$var1.pred


    y_hat_mint <-
        regression_kriging(
            daily_temperature_data$train,
            daily_temperature_data$test,
            daily_temperature_data$train[c("east", "north")],
            daily_temperature_data$test[c("east", "north")],

            regression_method = Random_Forest,
            formula = f_mint,

            cutoff = cutoff,
            width = width,
            vgm_model = vgm_model,

            nmax = Inf,
            flex_fit = TRUE,
            debug.level = 0,

            mtry = 5,
            min.node.size = 1
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
        "Regression Kriging Exp Random Forest",
        temperature_data$test$maxt,
        temperature_data$test$predicted_maxt,
        temperature_data$test$mint,
        temperature_data$test$predicted_mint
    )

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "RK_RF.csv"
#     ),
#     row.names = FALSE
# )
