
## Validation tests for daily rainfall grids ##
## Using the 80/20 split, I've established a collection
## of the best methods. Now we'll do a more thorough
## 10-fold CV to validate them fully

## At the moment, we have
## IDW - no regression
## Regression Kriging - using MLR
## Likelihood-based Kriging - with a trend

# %%

## Load starting data
load("Data/Daily_Rainfall/train_test_ten_fold_2016-2025.RData")

# %%

## Set up a column matrix, where we will save all
## the results for each method

library(dplyr)

final_results <- rain_data |>
    dplyr::select(stno, year, month, day, rain, fold)

## Example test - daily mean
daily_means <- rain_data |>
    group_by(year, month, day) |>
    mutate(predicted_rain = mean(rain))

final_results <- bind_cols(
    final_results,
    data.frame(`Daily Mean` = daily_means$predicted_rain)
)

# %%

## IDW

source("interpolation/IDW.R")

for (fold_idx in seq_along(unique(rain_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        daily_rain_data <-
            get_daily_rain_data(
                rain_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat <-
            IDW(
                daily_rain_data$train$y,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],
                idp = 2,
                nmax = 15
            )$var1.pred


        rain_data <-
            update_daily_predictions(
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
        "IDW",
        rain_data
    )

print(metrics_table)

write.csv(
    metrics_table,
    "Results/Daily_Rainfall/ten_fold_CV_2016-2025/IDW.csv",
    row.names = FALSE
)

final_results <- bind_cols(
    final_results,
    data.frame(`IDW` = rain_data$predicted_rain)
)

# %%

## Regression Kriging

source("regression/MLR.R")
source("interpolation/regression_kriging.R")

f <- y ~ east + north

for (fold_idx in seq_along(unique(rain_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        daily_rain_data <-
            get_daily_rain_data(
                rain_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat <-
            regression_kriging(
                daily_rain_data$train,
                daily_rain_data$test,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],

                regression_method = MLR,
                formula = f,
                cutoff = 350000,
                width = 15000,
                nmax = 15,

                flex_fit = TRUE,
                vgm_model = "Exp",
                debug.level = 0
            )$var1.pred


        rain_data <-
            update_daily_predictions(
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
    "Results/Daily_Rainfall/ten_fold_CV_2016-2025/RK_Exp_MLR.csv",
    row.names = FALSE
)

final_results <- bind_cols(
    final_results,
    data.frame(`Regression Kriging` = rain_data$predicted_rain)
)

# %%

## Likelihood-based Kriging (MLE)

source("interpolation/MLE_kriging.R")

f <- y ~ east + north

for (fold_idx in seq_along(unique(rain_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        daily_rain_data <-
            get_daily_rain_data(
                rain_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat <-
            MLE_kriging(
                daily_rain_data$train,
                daily_rain_data$test,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],


                cov_function = "Exp",

                formula = f,

                init_pars = c(0.001, 0.0001, 150),
                lower = c(1e-8, 0, 10),
                upper = c(0.05, 0.01, 500)
            )$pred

        rain_data <-
            update_daily_predictions(
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

write.csv(
    metrics_table,
    "Results/Daily_Rainfall/ten_fold_CV_2016-2025/MLE_Kriging.csv",
    row.names = FALSE
)

final_results <- bind_cols(
    final_results,
    data.frame(`MLE Kriging` = rain_data$predicted_rain)
)


# %%

f <- y ~ east + north + points5

for (fold_idx in seq_along(unique(rain_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        daily_rain_data <-
            get_daily_rain_data(
                rain_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat <-
            MLE_kriging(
                daily_rain_data$train,
                daily_rain_data$test,
                daily_rain_data$train[c("east", "north")],
                daily_rain_data$test[c("east", "north")],


                cov_function = "Exp",

                formula = f,

                init_pars = c(0.001, 0.0001, 150),
                lower = c(1e-8, 0, 10),
                upper = c(0.05, 0.01, 500)
            )$pred

        rain_data <-
            update_daily_predictions(
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

write.csv(
    metrics_table,
    "Results/Daily_Rainfall/ten_fold_CV_2016-2025/MLE_Kriging_points5.csv",
    row.names = FALSE
)

final_results <- bind_cols(
    final_results,
    data.frame(`MLE Kriging points5` = rain_data$predicted_rain)
)

# %%

# write.csv(
#     final_results,
#     "Results/Daily_Rainfall/ten_fold_CV_2016-2025/final_results.csv",
#     row.names = FALSE
# )
