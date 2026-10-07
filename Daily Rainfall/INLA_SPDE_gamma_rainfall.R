
## Validation tests for daily rainfall grids  ##

## INLA SPDE approach assuming the rainfall   ##
## data has a gamma distribution              ##

## There's a very ad-hoc thing in here where  ##
## I add +1 to y, then remove it from y_hat   ##
## This is to handle zeros, but isn't ideal   ##
## In future I'll use a hurdle model instead. ##

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
f <- y ~ 0 + Intercept +
    east + north +
    f(spatial.field, model = spde)

## Set up mesh
coast_points <- coordinates(island_outline)

ireland_nonconvex_hull <- fm_nonconvex_hull(
    coast_points,
    convex = -0.1,
)
ireland_mesh <- fm_mesh_2d_inla(
    boundary = ireland_nonconvex_hull,
    max.edge = c(10000, 40000),
    offset = c(10000, 40000),
    cutoff = 10000
)
ireland.spde <- inla.spde2.pcmatern(
    mesh = ireland_mesh,
    alpha = 1.5,
    prior.range = c(200000, 0.5),
    prior.sigma = c(0.04, 0.01)
)

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
    daily_rain_data$train$y <- daily_rain_data$train$y + 1

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
            data.frame(
                Intercept = 1,
                east = daily_rain_data$train$east,
                north = daily_rain_data$train$north
            )
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
            data.frame(
                Intercept = 1,
                east = daily_rain_data$test$east,
                north = daily_rain_data$test$north
            )
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
        family = "Gamma",
        control.predictor = list(
            A = inla.stack.A(ireland.stack),
            compute = FALSE,
            link = 1
        ),
        quantiles = NULL,
        control.compute = list(
            return.marginals = FALSE
        )
    )

    index.pred <- inla.stack.index(ireland.stack, "ireland.test")$data

    y_hat <- model$summary.fitted.values[index.pred, "mean"]
    y_hat <- y_hat - 1

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "INLA SPDE Gamma",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "INLA_SPDE_Gamma.csv"
#     ),
#     row.names = FALSE
# )

print(paste("Final Time", Sys.time() - st))


# %%

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
    daily_rain_data$y <- daily_rain_data$y + 1
    pred_grid <- grid_geodata

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
            data.frame(
                Intercept = 1,
                east = daily_rain_data$east,
                north = daily_rain_data$north
            )
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
            data.frame(
                Intercept = 1,
                east = pred_grid$east,
                north = pred_grid$north
            )
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
    y_hat <- y_hat - 1
    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        as.data.frame(daily_rain_data),
        plot_destination = paste0(
            "Figures/Daily_Rainfall/",
            "INLA_SPDE_gamma/INLA_SPDE_gamma_",
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