
## Validation tests for daily rainfall grids  ##

## INLA SPDE approach ##

## Ultimately this is still kriging             ##
## but the spatial model is with using Bayesian ##
## inference using the INLA-SPDE approach       ##

# %%

## Load in INLA library (I want to make a gridding function for this)
library(INLA)
library(fmesher)
library(dplyr)
library(sp)

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## Set formula
f <- "y ~ 0 + Intercept + f(spatial.field, model = spde) +"
linear_terms <- paste(
    "east", "north",
    "points5",
    # "dist2c",
    # "exp25k",
    # "n5", "e5", "s5", "w5",

    sep = " + "
)
f <- as.formula(paste(f, linear_terms))
linear_terms <- strsplit(linear_terms, " + ", fixed = TRUE)[[1]]

# %%

## Try setting up the mesh first
## since it's the same every time

coast_points <- coordinates(island_outline)

# ireland_convex_hull <- chull(coast_points)
# ireland_convex_hull <- coast_points[
#     c(ireland_convex_hull, ireland_convex_hull[1]),
# ]

ireland_nonconvex_hull <- fm_nonconvex_hull(
    coast_points,
    convex = -0.1,
)
ireland_mesh <- fm_mesh_2d_inla(
    boundary = ireland_nonconvex_hull,
    max.edge = c(20000, 60000),
    offset = c(20000, 60000),
    cutoff = 18000
)
ireland.spde <- inla.spde2.matern(mesh = ireland_mesh, alpha = 1.5)

s.index <- inla.spde.make.index(
    name = "spatial.field",
    n.spde = ireland.spde$n.spde
)

# %%

st <- Sys.time()

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    ## Get fixed effects from formula (covariates)
    fixed_effects <- lapply(
        linear_terms,
        function(x) daily_rain_data$train[[x]]
    )
    fixed_effects_test <- lapply(
        linear_terms,
        function(x) daily_rain_data$test[[x]]
    )
    names(fixed_effects_test) <- linear_terms
    names(fixed_effects) <- linear_terms
    fixed_effects$Intercept <- 1
    fixed_effects_test$Intercept <- 1

    coordinates(daily_rain_data$train) <- c("east", "north")
    coordinates(daily_rain_data$test) <- c("east", "north")
    proj4string(daily_rain_data$train) <- CRS("EPSG:29903")
    proj4string(daily_rain_data$test) <- CRS("EPSG:29903")

    A.train <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(daily_rain_data$train)
    )

    ## Stack for training data
    ireland.train.stack <- inla.stack(
        data = list(y = daily_rain_data$train$y),
        A = list(A.train, 1),
        effects = list(
            s.index,
            fixed_effects
        ),
        tag = "ireland.train"
    )

    ## Stack for test data
    A.test <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(daily_rain_data$test)
    )
    ireland.test.stack <- inla.stack(
        data = list(y = NA),
        A = list(A.test, 1),
        effects = list(
            s.index,
            fixed_effects_test
        ),
        tag = "ireland.test"
    )

    ## Full stack
    ireland.stack <- inla.stack(
        ireland.train.stack,
        ireland.test.stack
    )

    model <- inla(
        f,
        data = inla.stack.data(ireland.stack, spde = ireland.spde),
        family = "gaussian",
        control.predictor = list(
            A = inla.stack.A(ireland.stack),
            compute = FALSE
        ),
        quantiles = NULL,
        control.compute = list(
            return.marginals = FALSE
        )
    )

    index.pred <- inla.stack.index(ireland.stack, "ireland.test")$data

    y_hat <- model$summary.fitted.values[index.pred, "mean"]

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "INLA SPDE Exp",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "INLA_SPDE_Exp_trend_nonconvex.csv"
#     ),
#     row.names = FALSE
# )

print(paste("Final Time", Sys.time() - st))


# %%

## Check grids and computation time
times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    st <- Sys.time()

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index,
        experiment_type = "all_data"
    )
    pred_grid <- grid_geodata

    ## Get fixed effects from formula (covariates)
    fixed_effects <- lapply(
        linear_terms,
        function(x) daily_rain_data[[x]]
    )
    fixed_effects_test <- lapply(
        linear_terms,
        function(x) pred_grid[[x]]
    )
    names(fixed_effects_test) <- linear_terms
    names(fixed_effects) <- linear_terms
    fixed_effects$Intercept <- 1
    fixed_effects_test$Intercept <- 1

    coordinates(daily_rain_data) <- c("east", "north")
    coordinates(pred_grid) <- c("east", "north")
    proj4string(daily_rain_data) <- CRS("EPSG:29903")
    proj4string(pred_grid) <- CRS("EPSG:29903")

    A.train <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(daily_rain_data)
    )

    ## Stack for training data
    ireland.train.stack <- inla.stack(
        data = list(y = daily_rain_data$y),
        A = list(A.train, 1),
        effects = list(
            s.index,
            fixed_effects
        ),
        tag = "ireland.train"
    )

    ## Stack for test data
    A.test <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(pred_grid)
    )
    ireland.test.stack <- inla.stack(
        data = list(y = NA),
        A = list(A.test, 1),
        effects = list(
            s.index,
            fixed_effects_test
        ),
        tag = "ireland.test"
    )

    ## Full stack
    ireland.stack <- inla.stack(
        ireland.train.stack,
        ireland.test.stack
    )

    model <- inla(
        f,
        data = inla.stack.data(ireland.stack, spde = ireland.spde),
        family = "gaussian",
        control.predictor = list(
            A = inla.stack.A(ireland.stack),
            compute = FALSE
        ),
        quantiles = NULL,
        control.compute = list(
            return.marginals = FALSE
        )
    )

    index.pred <- inla.stack.index(ireland.stack, "ireland.test")$data

    y_hat <- model$summary.fitted.values[index.pred, "mean"]
    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        as.data.frame(daily_rain_data),
        plot_destination = paste0(
            "Figures/Daily_Rainfall/",
            "INLA_SPDE_Exp_trend/INLA_SPDE_Exp_trend_",
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
