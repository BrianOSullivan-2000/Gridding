## Validation tests for monthly rainfall grids ##
## Baseline regression kriging method we currently use ##

## So it's not QUITE the original baseline in this script
## because I'm using the 1991-2020 LTAs instead of
## 1981-2010. But given were testing from 2016-2025,
## I would expect the 1991-2020 LTAs to be better anyways

# %%

library(dplyr)

## Load starting data
load("Data/Monthly_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## Regression model
full_formula <- Nrr ~ east + I(east^2) + north + I(north^2) + east * north +
    n5 + s5 + e5 + w5 + ne5 + se5 + sw5 + nw5 + elev +
    exp5k + exp10k + exp15k + exp20k + exp25k + dist2c


# %%
## Gridding function

compute_monthly_grid <- function(train_data, target_data) {

    require(sp)
    require(geoR)
    require(gstat)

    train_data <- train_data[!duplicated(train_data[c("east", "north")]), ]

    current_month <- train_data[1, ]$month

    # Fit regression model for different exposures to the sea
    exp_models <- list()
    exp_models[[1]] <- lm(
        update(
            full_formula, . ~ . - exp10k - exp15k - exp20k - exp25k
        ), train_data
    )
    exp_models[[2]] <- lm(
        update(
            full_formula, . ~ . - exp5k - exp15k - exp20k - exp25k
        ), train_data
    )
    exp_models[[3]] <- lm(
        update(
            full_formula, . ~ . - exp5k - exp10k - exp20k - exp25k
        ), train_data
    )
    exp_models[[4]] <- lm(
        update(
            full_formula, . ~ . - exp5k - exp10k - exp15k - exp25k
        ), train_data
    )
    exp_models[[5]] <- lm(
        update(
            full_formula, . ~ . - exp5k - exp10k - exp15k - exp20k
        ), train_data
    )

    # Apply stepwise regression to each formula
    exp_models <- lapply(exp_models, function(model) {
        model$call$data <- quote(train_data)
        step(model, trace = 0)
    })

    # Pick model that returns highest adjusted R-squared value
    adjrsqs <- lapply(exp_models, {
        function(x) summary(x)$adj.r.squared
    })
    coeff <- as.data.frame(exp_models[[which.max(adjrsqs)]]$coefficients)
    names(coeff) <- "Beta"

    # Add residuals to data
    train_data$residuals <-
        as.numeric(exp_models[[which.max(adjrsqs)]]$residuals)

    # Identify outliers (defined as 2.5 * standard deviation(residuals))
    outlier_idx <-
        train_data$residuals < (-2.5 * sd(train_data$residuals, na.rm = TRUE)) |
        train_data$residuals > (2.5 * sd(train_data$residuals, na.rm = TRUE))

    # Convert to spatial data for producing spatial variogram
    train_data_geo <- as.geodata(
        train_data,
        coords.col = which(names(train_data) %in% c("east", "north")),
        data.col = which(names(train_data) == "residuals")
    )
    train_data_no_outliers <- as.geodata(
        as.data.frame(train_data_geo)[!outlier_idx, ]
    )

    # Get empirical variogram of residuals
    widths <- c(0, 4000, 8000, seq(16000, 80000, by = 8000))

    # Fit a covariance function to the empirical variogram
    vgm_no_outliers <- variog(
        train_data_no_outliers,
        uvec = widths, message = 0
    )
    fit_vgm_no_outliers <- variofit(vgm_no_outliers, message = 0)

    # Occasionally geoR will fit a sill or range of zero
    # then just fit a nugget model
    if ((fit_vgm_no_outliers$cov.pars[1] < 1e-16) ||
            (fit_vgm_no_outliers$cov.pars[2] < 1e-16)) {
        # Using the fitted model, define a variogram in the gstat package
        gstat_vgm <- vgm(
            model = "Nug",
            psill = fit_vgm_no_outliers$cov.pars[1] +
                fit_vgm_no_outliers$nugget,
            range = 0
        )
    } else {
        # Using the fitted model, define a variogram in the gstat package
        gstat_vgm <- vgm(
            model = "Exp",
            psill = fit_vgm_no_outliers$cov.pars[1],
            nugget = 0, range = fit_vgm_no_outliers$cov.pars[2]
        )
    }

    # Convert data to sp dataframe for gstat
    train_data_sp <- train_data
    coordinates(train_data_sp) <- c("east", "north")

    # Do the same for the target data (either test set or grid)
    coordinates(target_data) <- c("east", "north")

    # Interpolate the residuals using kriging (20 nearest neighbours)
    grid_vals <- krige(
        residuals ~ 1, train_data_sp, target_data, gstat_vgm,
        nmax = 20, debug.level = 0
    )$var1.pred
    # Add the regression trend back to the interpolated residuals
    grid_vals <-
        predict(exp_models[[which.max(adjrsqs)]], target_data) +
        grid_vals

    # Reverse the effect of dividing by 1981-2010 LTAs
    grid_vals <- grid_vals * (target_data[[paste0("m_", current_month)]])

    grid_vals
}

# %%

## Run CV

for (date_index in seq_len(nrow(dates))) {
    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    monthly_rain_data$train <- monthly_rain_data$train |>
        mutate(Nrr = rain / LTA)
    monthly_rain_data$test <- monthly_rain_data$test |>
        mutate(Nrr = rain)

    y_hat <-
        compute_monthly_grid(
            train_data = monthly_rain_data$train,
            target_data = monthly_rain_data$test
        )

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index,
        baseline = TRUE
    )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "Baseline Regression Kriging",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )

# write.csv(
#     metrics_table,
#     "Results/Monthly_Rainfall/train_test_80_20_2016-2025/baseline.csv",
#     row.names = FALSE
# )

# %%

## Plot grids and check computation time

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    st <- Sys.time()

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index,
        experiment_type = "all_data"
    ) |>
        mutate(Nrr = rain / LTA)

    target_data <- grid_geodata
    target_data[paste0("m_", current_month)] <-
        grid_LTAs_9120[paste0("m_", current_month)]

    y_hat <-
        compute_monthly_grid(
            train_data = monthly_rain_data,
            target_data = target_data
        )

    grid_geodata$rain <- y_hat

    monthly_rain_plot(
        grid_geodata,
        monthly_rain_data,
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/baseline/baseline_",
            current_year, "_",
            sprintf("%02d", current_month),
            ".jpg"
        )
    )

    et <- Sys.time()
    times <- c(times, et - st)
}

print(
    paste(
        "Mean Time",
        mean(times)
    )
)
