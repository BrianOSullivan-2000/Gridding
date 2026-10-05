
## Validation tests for monthly rainfall grids  ##

## INLA SPDE approach ##

## Ultimately this is still kriging             ##
## but the spatial model is with using Bayesian ##
## inference using the INLA-SPDE approach       ##

# %%

## Load in INLA library(
## I want to make a gridding function for this)
library(INLA)
library(fmesher)
library(dplyr)
library(sp)

## Load starting data
load("Data/Monthly_Rainfall/train_test_80_20_2016-2025.RData")

# %%

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

# %%

## Try setting up the mesh first
## since it's the same every time

coast_points <- coordinates(island_outline)

ireland_convex_hull <- chull(coast_points)
ireland_convex_hull <- coast_points[
    c(ireland_convex_hull, ireland_convex_hull[1]),
]
ireland_mesh <- fm_mesh_2d_inla(
    loc.domain = ireland_convex_hull,
    max.edge = c(40000, 80000),
    offset = c(40000, 80000),
    cutoff = 20000
)
ireland.spde <- inla.spde2.matern(mesh = ireland_mesh, alpha = 2)

s.index <- inla.spde.make.index(
    name = "spatial.field",
    n.spde = ireland.spde$n.spde
)

# %%

st <- Sys.time()

for (date_index in seq_len(nrow(dates))) {

    monthly_rain_data <- get_monthly_rain_data(
        rain_data, dates, date_index
    )

    ## Get fixed effects from formula (covariates)
    linear_terms <- all.vars(f[[3]])
    fixed_effects <- lapply(
        linear_terms,
        function(x) monthly_rain_data$train[[x]]
    )
    fixed_effects_test <- lapply(
        linear_terms,
        function(x) monthly_rain_data$test[[x]]
    )
    names(fixed_effects_test) <- linear_terms
    names(fixed_effects) <- linear_terms

    coordinates(monthly_rain_data$train) <- c("east", "north")
    coordinates(monthly_rain_data$test) <- c("east", "north")
    proj4string(monthly_rain_data$train) <- CRS("EPSG:29903")
    proj4string(monthly_rain_data$test) <- CRS("EPSG:29903")

    A.train <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(monthly_rain_data$train)
    )

    ## Stack for training data
    ireland.train.stack <- inla.stack(
        data = list(y = monthly_rain_data$train$y),
        A = list(A.train, 1),
        effects = list(
            c(s.index, list(Intercept = 1)),
            fixed_effects
        ),
        tag = "ireland.train"
    )

    ## Stack for test data
    A.test <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(monthly_rain_data$test)
    )
    ireland.test.stack <- inla.stack(
        data = list(y = NA),
        A = list(A.test, 1),
        effects = list(
            c(s.index, list(Intercept = 1)),
            fixed_effects_test
        ),
        tag = "ireland.test"
    )

    ## Full stack
    ireland.stack <- inla.stack(
        ireland.train.stack,
        ireland.test.stack
    )

    f_spatial <- update(
        f,
        . ~ . + Intercept - 1 + f(spatial.field, model = spde)
    )

    model <- inla(
        f_spatial,
        data = inla.stack.data(ireland.stack, spde = ireland.spde),
        family = "gaussian",
        control.predictor = list(
            A = inla.stack.A(ireland.stack), compute = TRUE
        ),
        control.compute = list(
            cpo = TRUE, dic = TRUE
        )
    )

    index.pred <- inla.stack.index(ireland.stack, "ireland.test")$data

    y_hat <- model$summary.fitted.values[index.pred, "mean"]

    rain_data <- update_monthly_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(date_index)
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "INLA SPDE",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Monthly_Rainfall/train_test_80_20_2016-2025/",
#         "INLA_SPDE_trend.csv"
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
    pred_grid <- grid_geodata

    coordinates(monthly_rain_data) <- c("east", "north")
    coordinates(pred_grid) <- c("east", "north")
    proj4string(monthly_rain_data) <- CRS("EPSG:29903")
    proj4string(pred_grid) <- CRS("EPSG:29903")

    A.train <- inla.spde.make.A(
        mesh = ireland_mesh,
        loc = coordinates(monthly_rain_data)
    )

    ## Stack for training data
    ireland.train.stack <- inla.stack(
        data = list(y = monthly_rain_data$y),
        A = list(A.train),
        effects = list(
            c(s.index, list(Intercept = 1))
            ## I have to add like the trend terms here
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
        A = list(A.test),
        effects = list(
            c(s.index, list(Intercept = 1))
        ),
        tag = "ireland.test"
    )

    ## Full stack
    ireland.stack <- inla.stack(
        ireland.train.stack,
        ireland.test.stack
    )

    f_temp <- y ~ -1 + Intercept + f(spatial.field, model = spde)

    model <- inla(
        f_temp,
        data = inla.stack.data(ireland.stack, spde = ireland.spde),
        family = "gaussian",
        control.predictor = list(
            A = inla.stack.A(ireland.stack), compute = TRUE
        ),
        control.compute = list(
            cpo = TRUE, dic = TRUE
        )
    )

    index.pred <- inla.stack.index(ireland.stack, "ireland.test")$data

    y_hat <- model$summary.fitted.values[index.pred, "mean"]
    grid_geodata$rain <-
        (y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    monthly_rain_plot(
        grid_geodata,
        as.data.frame(monthly_rain_data),
        plot_destination = paste0(
            "Figures/Monthly_Rainfall/INLA_SPDE/INLA_SPDE_",
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
