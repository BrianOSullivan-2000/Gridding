
## Validation tests for daily temperature grids ##
## Inverse Distance Weighting ##
## Simple method that uses weighted average based on inverse distance ##

## For spatio-temporal data, consider time as a third dimension, scaled
## using the hhyperparameter C

# %%

library(dplyr)
source("interpolation/spatio_temporal_kriging.R")
source("interpolation/ordinary_kriging.R")
source("utilities/variogram.R")

## Load starting data
load("Data/Daily_Temperature/train_test_five_fold_2025.RData")

# %%

## First just spatial kriging
for (fold_idx in seq_along(unique(temperature_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    for (date_index in seq_len(nrow(dates))) {

        daily_temperature_data <-
            get_daily_temperature_data(
                temperature_data,
                dates,
                date_index,
                current_fold = fold_idx
            )

        y_hat_maxt <-
            ordinary_kriging(
                daily_temperature_data$train$y_maxt,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                vgm_model = "Sph",
                cutoff = 400000,
                width = 20000,
                nmax = Inf,
                debug.level = 0
            )$var1.pred

        y_hat_mint <-
            ordinary_kriging(
                daily_temperature_data$train$y_mint,
                daily_temperature_data$train[c("east", "north")],
                daily_temperature_data$test[c("east", "north")],
                vgm_model = "Sph",
                cutoff = 400000,
                width = 20000,
                nmax = Inf,
                debug.level = 0
            )$var1.pred

        temperature_data <-
            update_daily_predictions(
                temperature_data,
                y_hat_maxt,
                y_hat_mint,
                dates,
                date_index,
                current_fold = fold_idx
            )
    }

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "Ordinary Kriging Sph",
        temperature_data |>
            as.data.frame()
    )

print(metrics_table)

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/OK_Sph.csv",
#     row.names = FALSE
# )

# %%

vgm_model <- "Sph"
cutoff <- 400000
width <- 20000
nmax <- Inf

plot_grid <- grid_geodata

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    ## First spatial IDW
    daily_temperature_data <-
        get_daily_temperature_data(
            temperature_data,
            dates,
            date_index,
            experiment_type = "all_data"
        )

    y_hat_maxt <-
        ordinary_kriging(
            daily_temperature_data$y_maxt,
            daily_temperature_data[c("east", "north")],
            plot_grid[c("east", "north")],
            vgm_model = vgm_model,
            cutoff = cutoff,
            width = width,
            nmax = nmax
        )$var1.pred

    y_hat_mint <-
        ordinary_kriging(
            daily_temperature_data$y_mint,
            daily_temperature_data[c("east", "north")],
            plot_grid[c("east", "north")],
            vgm_model = vgm_model,
            cutoff = cutoff,
            width = width,
            nmax = nmax
        )$var1.pred


    plot_grid$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    plot_grid$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        plot_grid,
        daily_temperature_data,
        plot_destination = paste0(
            "Figures/Daily_Temperature/five_fold_tests/",
            "OK_Sph/OK_Sph_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )
}

# %%

## Now spatio-temporal kriging
## Setting it all up

## Model parameters
nmax <- Inf
s_cutoff <- 400000
s_width <- 20000
t_cutoff <- 7
t_width <- 1

vgm_model <- "sumMetric"

spatial_vgm_model <- "Sph"
s_pars <- c(1, 500000, 0.1)
temporal_vgm_model <- "Mat"
t_pars <- c(20, 4, 0, 1)
joint_vgm_model <- "Exp"
j_pars <- c(1, 200000, 0, 0.5)

stAni <- 10000

parscales <- c(
    0.5, 10000, 0.1,
    0.5, 1, 0.1,
    0.5, 10000, 0.1,
    5000
)
lower_limits <- c(
    0.1, 5000, 0,
    0.1, 0.7, 0,
    0.1, 5000, 0,
    5000
)
upper_limits <- c(
    10, 500000, 1,
    30, 20, 1,
    10, 500000, 1,
    500000
)

maxit <- 400
vgm_progress <- TRUE
plot_vgm <- TRUE
return_emp_vgm <- FALSE
no_fit <- FALSE

fit.method <- 6

daily_window <- 7

## Prediction function - To save lines of code!
predict_temperature <- function(train_data, test_data, ...) {

    responses <- c(
        y_hat_maxt = "y_maxt",
        y_hat_mint = "y_mint"
    )

    predictions <- lapply(names(responses), function(prediction_name) {
        spatio_temporal_kriging(
            train_data,
            test_data,
            response = responses[[prediction_name]],

            nmax = nmax,
            s_cutoff = s_cutoff, s_width = s_width,
            t_cutoff = t_cutoff, t_width = t_width,

            vgm_model = vgm_model,

            spatial_vgm_model = spatial_vgm_model,
            temporal_vgm_model = temporal_vgm_model,
            joint_vgm_model = joint_vgm_model,
            s_pars = s_pars, t_pars = t_pars, j_pars = j_pars,

            parscales = parscales,
            lower_limits = lower_limits, upper_limits = upper_limits,
            maxit = maxit,

            flex_fit = TRUE,

            vgm_progress = vgm_progress, plot_vgm = plot_vgm,
            return_emp_vgm = return_emp_vgm, no_fit = no_fit,

            fit.method = fit.method,
            ...
        )$var1.pred
    })

    names(predictions) <- names(responses)
    predictions
}

# %%

## Validation

for (fold_idx in seq_along(unique(temperature_data$fold))) {

    print(paste0("Starting Fold ", fold_idx))
    fold_time <- Sys.time()

    prepared_temperature_data <-
        get_temperature_data(
            temperature_data,
            current_fold = fold_idx
        )

    y_hat <- predict_temperature(
        prepared_temperature_data$train,
        prepared_temperature_data$test
    )
    y_hat_maxt <- y_hat$y_maxt
    y_hat_mint <- y_hat$y_mint

    temperature_data <-
        update_predictions(
            temperature_data,
            y_hat_maxt,
            y_hat_mint,
            current_fold = fold_idx
        )

    print(paste0(
        "Fold ", fold_idx, " completed in ",
        Sys.time() - fold_time
    ))
}

# write.csv(
#     metrics_table,
#     "Results/Daily_Temperature/train_test_five_fold_2025/ST_OK_separable.csv",
#     row.names = FALSE
# )

# %%

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    date_indices <- (date_index - daily_window):(date_index + daily_window)
    date_indices <- date_indices[
        date_indices >= 1 & date_indices <= 365
    ]

    prepared_temperature_data <-
        get_temperature_data(
            temperature_data,
            dates = dates,
            date_indices = date_indices,
            experiment_type = "all_data"
        )

    plot_grid <- grid_geodata
    plot_grid$t <- date_index
    plot_grid <- bind_rows(
        plot_grid,
        plot_grid[1, ]
    )
    plot_grid[nrow(plot_grid), ]$t <- date_index + 1

    y_hat <- predict_temperature(
        prepared_temperature_data,
        plot_grid
    )
    y_hat_maxt <- y_hat$y_hat_maxt
    y_hat_mint <- y_hat$y_hat_mint

    y_hat_maxt <- y_hat_maxt[plot_grid$t == date_index]
    y_hat_mint <- y_hat_mint[plot_grid$t == date_index]

    plot_grid <- plot_grid[plot_grid$t == date_index, ]
    plot_grid$maxt <-
        y_hat_maxt + LTAs_max_9120[[paste0("m_", current_month)]]

    plot_grid$mint <-
        y_hat_mint + LTAs_min_9120[[paste0("m_", current_month)]]

    daily_temperature_plot(
        plot_grid,
        prepared_temperature_data |>
            filter(
                day == current_day,
                month == current_month,
                year == current_year
            ),
        plot_destination = paste0(
            "Figures/Daily_Temperature/five_fold_tests/",
            "ST_OK_separable/ST_OK_separable_fit_well_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
            ".jpg"
        )
    )
}
