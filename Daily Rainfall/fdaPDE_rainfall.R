
## Validation tests for daily rainfall grids ##

## functional data analysis with partial differential ##
## equations (fdaPDE). Models rainfall using basis    ##
## functions (as well as a few covariates)

## https://doi.org/10.1111/insr.12444

## Including the PDE hasn't shown much success, so       ##
## using the "noPDE" version, which basically just       ##
## uses the standard Laplacian to smooth basis functions ##


# %%

library(dplyr)

source("interpolation/fdaPDE_noPDE.R")

## Load starting data
load("Data/Daily_Rainfall/train_test_80_20_2016-2025.RData")

# %%

## For fdaPDE, you need to set up a mesh. I'll do this at
## the beginning to avoid computation overhead. If we use this
## method operationally, I imagine this as a piece of metadata
## you would load in

## Mesh-making libraries
library(fdaPDE2)
library(RTriangle)
library(sf)
library(concaveman)

island_points <- st_as_sf(island_outline, coords = c("x", "y"), crs = 29903)
island_points <- st_transform(island_points, crs = 4326)

boundary <- concaveman(island_points, concavity = 10)
boundary <- st_buffer(boundary, dist = 5000)
boundary <- st_simplify(boundary, dTolerance = 5000)

boundary_nodes <- st_cast(x = boundary, "POINT", crs = 4326)
boundary_nodes <- st_coordinates(x = boundary_nodes)
boundary_nodes <- data.frame(
    lon = boundary_nodes[, 1],
    lat = boundary_nodes[, 2]
)
boundary_nodes <- boundary_nodes[-nrow(boundary_nodes), ]
boundary_segments <- cbind(
    seq_len(nrow(boundary_nodes)),
    c(2:nrow(boundary_nodes), 1)
)

boundary_pslg <- pslg(P = boundary_nodes, S = boundary_segments)
mesh_regular <- triangulate(
    p = boundary_pslg, a = 0.002, q = 30
)
mesh_regular$H <- matrix(data = numeric(0), ncol = 2)
mesh <- triangulation(
    nodes = mesh_regular$P, cells = mesh_regular$T,
    boundary = mesh_regular$PB
)

# %%

## Regression Formula
f <- as.formula(
    paste(
        "y ~",
        paste(
            # "east", "north",
            "points5",
            "dist2c",
            sep = " + "
        )
    )
)

# %%

for (date_index in seq_len(nrow(dates))) {

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index
    )

    y_hat <-
        fdaPDE_noPDE(
            daily_rain_data$train,
            daily_rain_data$test,
            daily_rain_data$train[c("east", "north")],
            daily_rain_data$test[c("east", "north")],
            mesh = mesh, formula = f, lambda = 1e-5
        )$pred

    rain_data <- update_daily_predictions(
        rain_data, y_hat, dates, date_index
    )

    print(paste("Date", date_index))
}

metrics_table <-
    collect_metrics(
        metrics_table = NULL,
        "fdaPDE",
        rain_data$test$rain,
        rain_data$test$predicted_rain
    )
print(metrics_table)

# write.csv(
#     metrics_table,
#     paste0(
#         "Results/Daily_Rainfall/train_test_80_20_2016-2025/",
#         "fdaPDE.csv"
#     ),
#     row.names = FALSE
# )

# %%

## Plot grids and check computation time

times <- c()

plot_grid <- grid_geodata[c("east", "north")]

for (date_index in seq_len(nrow(dates))) {

    st <- Sys.time()

    current_day <- dates[date_index, ]$day
    current_month <- dates[date_index, ]$month
    current_year <- dates[date_index, ]$year

    daily_rain_data <- get_daily_rain_data(
        rain_data, dates, date_index, experiment_type = "all_data"
    )

    y_hat <-
        fdaPDE_noPDE(
            daily_rain_data,
            grid_geodata,
            daily_rain_data[c("east", "north")],
            grid_geodata[c("east", "north")],
            mesh = mesh, formula = f, lambda = 1e-5
        )$pred

    grid_geodata$rain <-
        expm1(y_hat) * grid_LTAs_9120[paste0("m_", current_month)]

    daily_rain_plot(
        grid_geodata,
        daily_rain_data,
        plot_destination = paste0(
            "Figures/Daily_Rainfall/fdaPDE/fdaPDE_",
            current_year, "_",
            sprintf("%02d", current_month), "_",
            sprintf("%02d", current_day),
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
