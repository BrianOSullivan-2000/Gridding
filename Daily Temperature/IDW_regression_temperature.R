
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

basic_terms <- c("east", "north")
safe_terms <- c("east", "north", "points5", "dist2c")
complex_terms <- c(
    "east", "north",
    "points5",
    "dist2c", "exp25k",
    "n5", "e5", "s5", "w5",
    "aspect_10km"
)
full_terms <- c(
    "east", "north",
    "points5",
    "dist2c", "exp5k", "exp10k", "exp15k",
    "exp20k", "exp25k",
    "n5", "e5", "s5", "w5",
    "ne5", "nw5", "se5", "sw5",
    "aspect_10km"
)
formulae_terms <- list(
    "Basic" = basic_terms,
    "Safe" = safe_terms,
    "Complex" = complex_terms,
    "Full" = full_terms
)

formulae_maxt <- lapply(
    formulae_terms,
    function(terms) {
        as.formula(
            paste(
                "y_maxt ~",
                paste(
                    terms,
                    collapse = " + "
                )
            )
        )
    }
)

formulae_mint <- lapply(
    formulae_terms,
    function(terms) {
        as.formula(
            paste(
                "y_mint ~",
                paste(
                    terms,
                    collapse = " + "
                )
            )
        )
    }
)

## Regression Formula
f_maxt <- as.formula(
    paste(
        "y_maxt ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            # "exp25k",
            "aspect_10km",
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
            # "exp25k",
            "aspect_10km",
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

# %%

# Random Forests Grid Search

source("regression/random_forest.R")

## Set up grid search for hyperparameters
idps <- c(1, 2)
nmaxs <- c(10, 15, 20)
mtrys <- c(3, 5, 8, 10)
min.node.sizes <- c(1, 2, 3)
formula_idx <- seq_along(formulae_maxt)
hyperparameters <- expand.grid(
    idp = idps,
    nmax = nmaxs,
    mtry = mtrys,
    min.node.size = min.node.sizes,
    formula_idx = formula_idx
)

for (hyperparameter_index in seq_len(nrow(hyperparameters))) {

    idp <- hyperparameters[hyperparameter_index, ]$idp
    nmax <- hyperparameters[hyperparameter_index, ]$nmax
    mtry <- hyperparameters[hyperparameter_index, ]$mtry
    min.node.size <- hyperparameters[hyperparameter_index, ]$min.node.size
    formula_idx <- hyperparameters[hyperparameter_index, ]$formula_idx

    if (
        mtry <
            length(
                attr(
                    terms(formulae_maxt[[formula_idx]]),
                    "term.labels"
                )
            )
    ) {

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
                    regression_method = Random_Forest,
                    formula = formulae_maxt[[formula_idx]],
                    mtry = mtry,
                    min.node.size = min.node.size,
                    idp = idp,
                    nmax = nmax
                )$var1.pred


            y_hat_mint <-
                regression_IDW(
                    daily_temperature_data$train,
                    daily_temperature_data$test,
                    daily_temperature_data$train[c("east", "north")],
                    daily_temperature_data$test[c("east", "north")],
                    regression_method = Random_Forest,
                    formula = formulae_mint[[formula_idx]],
                    mtry = mtry,
                    min.node.size = min.node.size,
                    idp = idp,
                    nmax = nmax
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
                    "IDW RF  ",
                    "idp-", idp,
                    " nmax-", nmax,
                    " mtry-", mtry,
                    " min.node.size-", min.node.size, " ",
                    c(
                        "basic", "safe", "complex", "full"
                    )[formula_idx]
                ),
                temperature_data$test$maxt,
                temperature_data$test$predicted_maxt,
                temperature_data$test$mint,
                temperature_data$test$predicted_mint
            )
    }

    print(hyperparameter_index)
}

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Temperature/train_test_80_20_2025/",
#         "IDW_RF_config.csv"
#     ),
#     row.names = FALSE
# )

## I've landed on idp = 1, nmax = 15, mtry = 5, min.node.size = 1,
## east + north + points5 + dist2c + exp25k +
## n5 + e5 + s5 + w5 + aspect_10km


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

    y_hat_maxt <-
        regression_IDW(
            daily_temperature_data,
            grid_geodata,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = MLR,
            formula = f_maxt,

            idp = 1,
            nmax = 15
        )$var1.pred

    y_hat_mint <-
        regression_IDW(
            daily_temperature_data,
            grid_geodata,
            daily_temperature_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            regression_method = MLR,
            formula = f_mint,

            idp = 1,
            nmax = 15
        )$var1.pred

    grid_geodata$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    grid_geodata$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        grid_geodata,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/IDW_MLR/",
            "IDW_MLR_",
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
