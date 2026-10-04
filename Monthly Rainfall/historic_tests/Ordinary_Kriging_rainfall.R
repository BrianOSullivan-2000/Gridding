
## Validation tests for monthly rainfall grids ##

## Ordinary Kriging ##
## Assumes data follows a Gaussian process with a spatial dependence ##
## Data is predicted according to distance from nearby points ##

# %%

library(dplyr)
source("interpolation/ordinary_kriging.R")

## Load starting data
load("Data/Monthly_Rainfall/train_test_80_20_1883-1940.RData")

# %%

## Also going to check against nearest neighbours
nmaxs <- c(8, 12, 15, 20, 50, Inf)

for (hyperparameter_index in seq_along(nmaxs)) {

    nmax <- nmaxs[hyperparameter_index]

    for (date_index in seq_len(nrow(dates))) {

        monthly_rain_data <- get_monthly_rain_data(
            rain_data, dates, date_index
        )

        y_hat <-
            ordinary_kriging(
                monthly_rain_data$train$y,
                monthly_rain_data$train[c("east", "north")],
                monthly_rain_data$test[c("east", "north")],

                cutoff = 350000,
                width = 15000,

                nmax = nmax,

                flex_fit = TRUE,
                vgm_model = "Mat",
                debug.level = 0
            )$var1.pred

        rain_data <- update_monthly_predictions(
            rain_data, y_hat, dates, date_index
        )
    }

    metrics_table <-
        collect_metrics(
            metrics_table,
            paste0("OK nmax-", nmax),
            rain_data$test$rain,
            rain_data$test$predicted_rain
        )

    # print(metrics_table)
    print(hyperparameter_index)
}

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Monthly_Rainfall/train_test_80_20_1883-1940/",
#         "OK_Mat_config_nmax.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Check grids and computation time

## Exp with nmax = 20 and cutoff = 350000

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    st <- Sys.time()

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        ordinary_kriging(
            monthly_rain_data$y,
            monthly_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],

            cutoff = 350000,
            width = 15000,
            nmax = 20,

            flex_fit = TRUE,
            vgm_model = "Exp",
            debug.level = 0
        )$var1.pred

    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    monthly_rain_plot(
        grid_geodata,
        monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/historic/",
            "OK_Exp/OK_Exp_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)

    print(date_index)
}

print("Mean Times")
print(mean(times))
