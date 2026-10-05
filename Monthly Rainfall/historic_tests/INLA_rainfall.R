
## Validation tests for monthly rainfall grids  ##

## INLA SPDE approach ##

## Ultimately this is still kriging             ##
## but the spatial model is with using Bayesian ##
## inference using the INLA-SPDE approach       ##

# %%

## Load in INLA library (I want to make a gridding function for this)
library(INLA)
library(dplyr)

## Load 2016-2025 data & get rainfall
load("Data/Monthly_Rainfall/train_test_80_20_2016-2025.RData")
rain_2016_2025 <- rain_data
dates_2016_2025 <- dates
gmrd_2016_2025 <- get_monthly_rain_data

## Load starting data
load("Data/Monthly_Rainfall/train_test_80_20_1883-1940.RData")

# %%

## Going to use INLA in the regression step

## Set formula
f <- as.formula(
    paste(
        "y ~",
        paste(
            "east", "north",
            "points5",
            "dist2c",
            # "exp25k",
            # "n5", "e5", "s5", "w5",
            sep = " + "
        )
    )
)

## Pick out priors
## I'll do this in an empirical way,
## Run a regression over every month (2016-2025) after normalising with
## LTAs. We can see how the betas look for setting priors

betas <-
    matrix(
        nrow = 0,
        ncol = length(attr(terms(f), "term.labels")) + 1
    )
colnames(betas) <- c(
    "(Intercept)",
    attr(terms(f), "term.labels")
)

for (date_index in seq_len(nrow(dates_2016_2025))) {
    monthly_rain_data <- gmrd_2016_2025(
        rain_2016_2025, dates_2016_2025, date_index
    )

    model <- lm(
        f, monthly_rain_data$train
    )
    betas <- rbind(
        betas,
        model$coefficients
    )
}

# %%

source("interpolation/ordinary_kriging.R")


st <- Sys.time()

for (date_index in seq_len(nrow(dates))) {

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    ## Train linear model using INLA
    model <- inla(
        formula = f,
        data = monthly_rain_data$train,

        ## Controls for priors
        ## i.e. what we think the trends are, and
        ## how strongly we believe them
        control.fixed = list(
            mean.intercept = mean(betas[, 1]),
            prec.intercept = 4,
            mean = as.list(apply(betas, 2, mean)),
            prec = as.list(
                apply(
                    betas, 2,
                    function(b) 10^floor(log10(abs(mean(b))) + 4)
                )[2:ncol(betas)]
            )
        )
    )

    ## Get model residuals
    monthly_betas <- model$summary.fixed[, 1]
    y_trend <-
        model.matrix(f, monthly_rain_data$train) %*% monthly_betas
    y_res <- monthly_rain_data$train$y - y_trend

    ## Interpolate y_hat
    y_hat <- ordinary_kriging(
        y_res,
        monthly_rain_data$train[c("east", "north")],
        monthly_rain_data$test[c("east", "north")],
        nmax = 20,
        cutoff = 350000,
        width = 15000,
        flex_fit = TRUE,
        vgm_model = "Exp",
        debug.level = 0
    )$var1.pred

    y_hat <- y_hat +
        model.matrix(f, monthly_rain_data$test) %*% monthly_betas

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "RK Exp INLA",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Monthly_Rainfall/train_test_80_20_1883-1940/",
#         "RK_Exp_INLA_strong_prior.csv"
#     ),
#     row.names = FALSE
# )

print(paste("Final Time", Sys.time() - st))

# %%

## Check grids and computation time
times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    st <- Sys.time()

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index,
        experiment_type = "all_data"
    )

    ## Train linear model using INLA
    model <- inla(
        formula = f,
        data = monthly_rain_data,

        ## Controls for priors
        ## i.e. what we think the trends are, and
        ## how strongly we believe them
        control.fixed = list(
            mean.intercept = mean(betas[, 1]),
            prec.intercept = 1,
            mean = as.list(apply(betas, 2, mean)),
            prec = as.list(
                apply(
                    betas, 2,
                    function(b) 10^floor(log10(abs(mean(b))))
                )[2:ncol(betas)]
            )
        )
    )

    ## Get model residuals
    monthly_betas <- model$summary.fixed[, 1]
    y_trend <-
        model.matrix(f, monthly_rain_data) %*% monthly_betas
    y_res <- monthly_rain_data$y - y_trend

    ## Interpolate y_hat
    y_hat <- ordinary_kriging(
        y_res,
        monthly_rain_data[c("east", "north")],
        grid_geodata[c("east", "north")],
        nmax = 20,
        cutoff = 350000,
        width = 15000,
        flex_fit = TRUE,
        vgm_model = "Exp",
        debug.level = 0
    )$var1.pred

    grid_geodata$y <- 0
    y_hat <- y_hat +
        model.matrix(f, grid_geodata) %*% monthly_betas

    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    monthly_rain_plot(
        grid_geodata,
        monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/historic/",
            "RK_Exp_INLA/RK_Exp_INLA_",
            current_year, "_",
            sprintf("%02d", current_month),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)

    print(date_index)
}

print("Mean Times")
print(mean(times))
