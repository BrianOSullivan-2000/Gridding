
## Validation tests for monthly rainfall grids ##
## Using the 80/20 split, I've established a collection
## of the best methods. Now we'll do a more thorough
## 10-fold CV to validate them fully

## At the moment, we have
## IDW - with GLS regression
## Regression Kriging - using GLS
## Likelihood-based Kriging - with a trend

# %%

## Load starting data
load("Data/Monthly_Rainfall/train_test_ten_fold_2016-2025.RData")

# %%

## Set up a column matrix, where we will save all
## the results for each method

library(dplyr)

final_results <- rain_data |>
    dplyr::select(stno, year, month, rain, fold)

## Example test - monthly mean
monthly_means <- rain_data |>
    group_by(year, month) |>
    mutate(predicted_rain = mean(rain))

final_results <- bind_cols(
    final_results,
    data.frame(`Monthly Mean` = monthly_means$predicted_rain)
)

# %%

## Regression IDW using GLS for regression

source("interpolation/regression_IDW.R")
source("regression/GLS.R")

f <- y ~ east + north + points5

for (fold_idx in seq_along(unique(rain_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        monthly_rain_data <-
            get_monthly_rain_data(
                rain_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat <-
            regression_IDW(
                monthly_rain_data$train,
                monthly_rain_data$test,
                monthly_rain_data$train[c("east", "north")],
                monthly_rain_data$test[c("east", "north")],

                regression_method = GLS,
                formula = f,
                flex_fit = TRUE, flex_vgm = TRUE,
                vgm_model = "Exp",

                idp = 1.5,
                nmax = 15
            )$var1.pred


        rain_data <-
            update_monthly_predictions(
                rain_data,
                y_hat,
                dates,
                date_index,
                current_fold = fold_idx
            )

        print(paste0(
            "Date index ", date_index, " completed"
        ))
    }

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}
dates[date_index, ]

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "IDW",
        rain_data
    )

print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Monthly_Rainfall/ten_fold_CV_2016-2025/IDW_GLS.csv",
#     row.names = FALSE
# )

final_results <- bind_cols(
    final_results,
    data.frame(`IDW GLS` = rain_data$predicted_rain)
)

# %%

## Regression Kriging using GLS

source("interpolation/regression_kriging.R")

f <- y ~ east + north + points5

for (fold_idx in seq_along(unique(rain_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        monthly_rain_data <-
            get_monthly_rain_data(
                rain_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat <- tryCatch(
            {
                regression_kriging(
                    monthly_rain_data$train,
                    monthly_rain_data$test,
                    monthly_rain_data$train[c("east", "north")],
                    monthly_rain_data$test[c("east", "north")],
                    regression_method = GLS,
                    formula = f,
                    cutoff = 350000,
                    width = 15000,
                    nmax = 20,
                    flex_fit = TRUE,
                    vgm_model = "Exp",
                    debug.level = 0
                )$var1.pred
            },
            error = function(e) {
                message("First attempt failed: ", conditionMessage(e))
                message("Retrying with width = 20000")

                regression_kriging(
                    monthly_rain_data$train,
                    monthly_rain_data$test,
                    monthly_rain_data$train[c("east", "north")],
                    monthly_rain_data$test[c("east", "north")],
                    regression_method = GLS,
                    formula = f,
                    cutoff = 350000,
                    width = 20000,
                    nmax = 20,
                    flex_fit = TRUE,
                    vgm_model = "Exp",
                    debug.level = 0
                )$var1.pred
            }
        )

        rain_data <-
            update_monthly_predictions(
                rain_data,
                y_hat,
                dates,
                date_index,
                current_fold = fold_idx
            )

        print(paste0(
            "Date index ", date_index, " completed"
        ))
    }

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "Regression Kriging",
        rain_data
    )

print(metrics_table)

write.csv(
    metrics_table,
    "Results/Monthly_Rainfall/ten_fold_CV_2016-2025/RK_Exp_MLR.csv",
    row.names = FALSE
)

final_results <- bind_cols(
    final_results,
    data.frame(`Regression Kriging GLS` = rain_data$predicted_rain)
)

# %%

## Likelihood-based Kriging (MLE)

source("interpolation/MLE_kriging.R")

f <- y ~ east + north + east^2 + north^2 + points5

for (fold_idx in seq_along(unique(rain_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        monthly_rain_data <-
            get_monthly_rain_data(
                rain_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat <-
            MLE_kriging(
                monthly_rain_data$train |>
                    mutate(east = east / 1000, north = north / 1000),
                monthly_rain_data$test |>
                    mutate(east = east / 1000, north = north / 1000),
                monthly_rain_data$train[c("east", "north")],
                monthly_rain_data$test[c("east", "north")],

                cov_function = "Exp",

                formula = f,

                init_pars = c(0.001, 0.0001, 150),
                lower = c(1e-8, 0, 10),
                upper = c(0.05, 0.01, 500)
            )$pred

        rain_data <-
            update_monthly_predictions(
                rain_data,
                y_hat,
                dates,
                date_index,
                current_fold = fold_idx
            )

        print(paste0(
            "Date index ", date_index, " completed"
        ))
    }

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}


metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "Regression Kriging MLE",
        rain_data
    )

print(metrics_table)

write.csv(
    metrics_table,
    "Results/Monthly_Rainfall/ten_fold_CV_2016-2025/MLE_Kriging.csv",
    row.names = FALSE
)

final_results <- bind_cols(
    final_results,
    data.frame(`MLE Kriging` = rain_data$predicted_rain)
)
