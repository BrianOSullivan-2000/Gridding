
## This script is for investigating the covariance model fitting
## to spatio-temporal empirical variograms
## The main application here is for spatio-temporal kriging

# %%

library(dplyr)
source("utilities/variogram.R")

## Load starting data
load("Data/Daily_Temperature/train_test_five_fold_2025.RData")

# %%

## The spatio-temporal variogram has
## lots of fine-tuning and arguments needed

## Empirical variogram parameters
s_cutoff <- 400000
s_width <- 20000
t_cutoff <- 7
t_width <- 1

vgm_model <- "sumMetric"

spatial_vgm_model <- "Sph"
s_pars <- c(1, 500000, 0.1)
temporal_vgm_model <- "Mat"
t_pars <- c(50, 7, 0, 1.5)
joint_vgm_model <- "Exp"
j_pars <- c(1, 200000, 0, 0.5)

stAni <- 300000

parscales <- c(
    0.5, 10000, 0.1,
    0.5, 1, 0.1,
    0.5, 10000, 0.1,
    5000
)
lower_limits <- c(
    0.01, 1000, 0,
    0.01, 0.3, 0,
    0.01, 1000, 0,
    100
)
upper_limits <- c(
    50, 2000000, 5,
    200, 500, 5,
    50, 2000000, 5,
    5000000
)

maxit <- 1000
vgm_progress <- TRUE
plot_vgm <- TRUE
return_emp_vgm <- FALSE
no_fit <- FALSE

fit.method <- 10
daily_window <- 7

fit_ST_variograms <- function(df, ...) {

    vgmsts <- list()

    for (response in c("y_maxt", "y_mint")) {

        df$y <- df[[response]]

        vgmsts[[response]] <-
            spatio_temporal_variogram(
                df,

                s_cutoff = s_cutoff, s_width = s_width,
                t_cutoff = t_cutoff, t_width = t_width,

                vgm_model = vgm_model,

                spatial_vgm_model = spatial_vgm_model,
                temporal_vgm_model = temporal_vgm_model,
                joint_vgm_model = joint_vgm_model,
                # s_pars = s_pars, t_pars = t_pars, j_pars = j_pars,

                parscales = parscales,
                lower_limits = lower_limits, upper_limits = upper_limits,
                maxit = maxit,

                flex_fit = TRUE,

                vgm_progress = vgm_progress, plot_vgm = plot_vgm,
                return_emp_vgm = return_emp_vgm, no_fit = no_fit,

                fit.method = fit.method,
                ...
            )$vgmst

    }

    vgmsts
}

# %%

# for (date_index in seq_len(nrow(dates))) {
for (date_index in seq(10, 300, by = 20)) {

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

    vgms <-
        fit_ST_variograms(
            prepared_temperature_data
        )
}

vgms
