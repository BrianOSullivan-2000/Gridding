
## Validation tests for daily rainfall grids  ##

## INLA SPDE approach using a hurdle model    ##
## So we're modelling occurrence of rain (i.e. ##
## a Bernoulli distribution), and then the    ##
## quantity of rain (a Gamma distribution)    ##

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
f <- Z ~ -1 + occurrence.intercept + quantity.intercept +
    f(
        occurrence.spatial.field,
        model = ireland.spde
    ) +
    f(
        occurrence.quantity.spatial.field,
        copy = "occurrence.spatial.field",
        fixed = FALSE,
        hyper = list(
            theta = list(prior = "gaussian", param = c(0, 1))
        )
    ) +
    f(
        quantity.spatial.field,
        model = ireland.spde
    )

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
    prior.sigma = c(0.9, 0.01)
)

occurrence.index <- inla.spde.make.index(
    name = "occurrence.spatial.field",
    n.spde = ireland.spde$n.spde
)
occurrence.quantity.index <- inla.spde.make.index(
    name = "occurrence.quantity.spatial.field",
    n.spde = ireland.spde$n.spde
)
quantity.index <- inla.spde.make.index(
    name = "quantity.spatial.field",
    n.spde = ireland.spde$n.spde
)

# %%

st <- Sys.time()

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    ## Define occurrence variable
    daily_rain_data$train <- daily_rain_data$train |>
        mutate(
            occurrence = as.numeric(normalized_rain > 0),
            y = if_else(occurrence == 0, NA_real_, normalized_rain)
        )
    daily_rain_data$test <- daily_rain_data$test |>
        mutate(
            occurrence = as.numeric(normalized_rain > 0),
            y = if_else(occurrence == 0, NA_real_, normalized_rain)
        )

    coordinates(daily_rain_data$train) <- c("east", "north")
    coordinates(daily_rain_data$test) <- c("east", "north")
    proj4string(daily_rain_data$train) <- CRS("EPSG:29903")
    proj4string(daily_rain_data$test) <- CRS("EPSG:29903")

    A.train <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(daily_rain_data$train)
    )

    ## Stacks for training data
    ireland.train.occurrence.stack <- inla.stack(
        data = list(
            Z = cbind(
                as.vector(daily_rain_data$train$occurrence),
                NA
            ),
            link = 1
        ),
        A = list(A.train, 1),
        effects = list(
            occurrence.index,
            list(
                occurrence.intercept =
                    rep(1, nrow(daily_rain_data$train))
            )
        ),
        tag = "ireland.train.occurrence"
    )
    ireland.train.quantity.stack <- inla.stack(
        data = list(
            Z = cbind(
                NA,
                as.vector(daily_rain_data$train$y)
            ),
            link = 2
        ),
        A = list(A.train, 1),
        effects = list(
            c(occurrence.quantity.index, quantity.index),
            quantity.intercept =
                rep(1, nrow(daily_rain_data$train))
        ),
        tag = "ireland.train.quantity"
    )

    ## Stacks for test data
    A.test <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(daily_rain_data$test)
    )

    ireland.test.occurrence.stack <- inla.stack(
        data = list(
            Z = matrix(NA, nrow(A.test), 2),
            link = 1
        ),
        A = list(A.test, 1),
        effects = list(
            occurrence.index,
            occurrence.intercept =
                rep(1, nrow(daily_rain_data$test))
        ),
        tag = "ireland.test.occurrence"
    )
    ireland.test.quantity.stack <- inla.stack(
        data = list(
            Z = matrix(NA, nrow(A.test), 2),
            link = 2
        ),
        A = list(A.test, 1),
        effects = list(
            c(occurrence.quantity.index, quantity.index),
            quantity.intercept =
                rep(1, nrow(daily_rain_data$test))
        ),
        tag = "ireland.test.quantity"
    )

    ## Full stack
    ireland.stack <- inla.stack(
        ireland.train.occurrence.stack,
        ireland.train.quantity.stack,
        ireland.test.occurrence.stack,
        ireland.test.quantity.stack
    )
    link <- inla.stack.data(ireland.stack)$link

    model <- inla(
        f,
        family = c("binomial", "gamma"),
        data = inla.stack.data(ireland.stack, spde = ireland.spde),
        control.predictor = list(
            A = inla.stack.A(ireland.stack),
            compute = TRUE,
            link = link
        ),
        quantiles = NULL,
        control.compute = list(
            return.marginals = FALSE
        ),
        control.inla = list(
            strategy = "adaptive",
            int.strategy = "eb"
        )
    )

    occurrence.index.pred <- inla.stack.index(
        ireland.stack,
        "ireland.test.occurrence"
    )$data
    quantity.index.pred <- inla.stack.index(
        ireland.stack,
        "ireland.test.quantity"
    )$data

    y_hat_occurrence <-
        model$summary.fitted.values[occurrence.index.pred, "mean"]
    y_hat_quantity <-
        model$summary.fitted.values[quantity.index.pred, "mean"]

    y_hat <- y_hat_occurrence * y_hat_quantity

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

write.csv(
    metrics_table,
    paste0(
        "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
        "INLA_SPDE_hurdle.csv"
    ),
    row.names = FALSE
)

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

    ## Define occurrence variable
    daily_rain_data <- daily_rain_data |>
        mutate(
            occurrence = as.numeric(normalized_rain > 0),
            y = if_else(occurrence == 0, NA_real_, normalized_rain)
        )

    pred_grid <- grid_geodata |>
        mutate(occurrence = NA_real_, y = NA_real_)

    coordinates(daily_rain_data) <- c("east", "north")
    coordinates(pred_grid) <- c("east", "north")
    proj4string(daily_rain_data) <- CRS("EPSG:29903")
    proj4string(pred_grid) <- CRS("EPSG:29903")


    ## Stack for training data
    A.train <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(daily_rain_data)
    )

    ireland.train.occurrence.stack <- inla.stack(
        data = list(
            Z = cbind(
                as.vector(daily_rain_data$occurrence),
                NA
            ),
            link = 1
        ),
        A = list(A.train, 1),
        effects = list(
            occurrence.index,
            list(
                occurrence.intercept =
                    rep(1, nrow(daily_rain_data))
            )
        ),
        tag = "ireland.train.occurrence"
    )
    ireland.train.quantity.stack <- inla.stack(
        data = list(
            Z = cbind(
                NA,
                as.vector(daily_rain_data$y)
            ),
            link = 2
        ),
        A = list(A.train, 1),
        effects = list(
            c(occurrence.quantity.index, quantity.index),
            quantity.intercept =
                rep(1, nrow(daily_rain_data))
        ),
        tag = "ireland.train.quantity"
    )



    ## Stack for test data
    A.test <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(pred_grid)
    )
    ireland.test.occurrence.stack <- inla.stack(
        data = list(
            Z = matrix(NA, nrow(A.test), 2),
            link = 1
        ),
        A = list(A.test, 1),
        effects = list(
            occurrence.index,
            occurrence.intercept =
                rep(1, nrow(pred_grid))
        ),
        tag = "ireland.test.occurrence"
    )
    ireland.test.quantity.stack <- inla.stack(
        data = list(
            Z = matrix(NA, nrow(A.test), 2),
            link = 2
        ),
        A = list(A.test, 1),
        effects = list(
            c(occurrence.quantity.index, quantity.index),
            quantity.intercept =
                rep(1, nrow(pred_grid))
        ),
        tag = "ireland.test.quantity"
    )


    ## Full stack
    ireland.stack <- inla.stack(
        ireland.train.occurrence.stack,
        ireland.train.quantity.stack,
        ireland.test.occurrence.stack,
        ireland.test.quantity.stack
    )
    link <- inla.stack.data(ireland.stack)$link


    model <- inla(
        f,
        family = c("binomial", "gamma"),
        data = inla.stack.data(ireland.stack, spde = ireland.spde),
        control.predictor = list(
            A = inla.stack.A(ireland.stack),
            compute = TRUE,
            link = link
        ),
        quantiles = NULL,
        control.compute = list(
            return.marginals = FALSE
        ),
        control.inla = list(
            strategy = "adaptive",
            int.strategy = "eb"
        )
    )


    occurrence.index.pred <- inla.stack.index(
        ireland.stack,
        "ireland.test.occurrence"
    )$data
    quantity.index.pred <- inla.stack.index(
        ireland.stack,
        "ireland.test.quantity"
    )$data

    y_hat_occurrence <-
        model$summary.fitted.values[occurrence.index.pred, "mean"]
    y_hat_quantity <-
        model$summary.fitted.values[quantity.index.pred, "mean"]
    y_hat <- y_hat_occurrence * y_hat_quantity

    grid_geodata$rain <-
        y_hat * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        as.data.frame(daily_rain_data),
        plot_destination = paste0(
            "Figures/Daily_Rainfall/",
            "INLA_SPDE_hurdle/INLA_SPDE_hurdle_",
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