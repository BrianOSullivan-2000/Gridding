
## I'm happy with MLE Kriging as a findal gridding method ##
## for the monthly rainfall grids (modern + historic)     ##

## I'll start by making 1883-1940 first

# %%

## First read in our rainfall data
rainfall_data <- read.csv(
    "Data/Monthly_Rainfall/monthly_rain_ROI_NI_1883-1940.csv"
)
## Remove any missing values
rainfall_data <- rainfall_data[!is.na(rainfall_data$rain), ]

## From this we're getting geodata and LTAs
load("Data/Monthly_Rainfall/train_test_80_20_1883-1940.RData")

prediction_grid <- grid_geodata |>
    dplyr::select(east, north, points5)

## Finally, source the MLE function
source("interpolation/MLE_kriging.R")

# %%

compute_monthly_grid <- function(
    historic_rain,
    current_year,
    current_month
) {

    require(dplyr)

    ## Get data for current month
    monthly_rain_data <- historic_rain |>
        filter(year == current_year, month == current_month)

    ## Add geodata and LTAs
    monthly_rain_data <- monthly_rain_data |>
        left_join(station_geodata, by = "stno") |>
        filter(!is.na(east), !is.na(north), !is.na(points5)) |>
        left_join(station_LTAs, by = "stno") |>
        mutate(LTA = .data[[paste0("m_", current_month)]]) |>
        dplyr::select(
            stno, year, month, rain, east, north, points5, LTA
        )

    ## Transform data
    monthly_rain_data <- monthly_rain_data |>
        mutate(normalized_rain = rain / LTA)

    ## Drop duplicate locations
    monthly_rain_data <- monthly_rain_data[
        !duplicated(monthly_rain_data[c("east", "north")]),
    ]

    ## Regresison model
    f <- normalized_rain ~ east + north + points5

    ## Actually gridding (Kriging algorithm using Max. Likelihood)
    rain_pred <-
        MLE_kriging(
            monthly_rain_data,
            prediction_grid,
            monthly_rain_data[c("east", "north")],
            prediction_grid[c("east", "north")],

            cov_function = "Exp",

            formula = f,

            init_pars = c(0.001, 0.0001, 150),
            lower = c(1e-8, 0, 10),
            upper = c(0.05, 0.01, 500)
        )$pred

    ## This is a quick sanity check for negative values
    if (sum(rain_pred < 0) > 0) {
        warning(paste(
            sum(rain_pred < 0), "negative grid points detected in",
            month.abb[current_month],
            current_year, "- set to zero"
        ))
        rain_pred[rain_pred < 0] <- 0
    }

    ## Reverse transformation
    prediction_grid$rain <-
        rain_pred * grid_LTAs[paste0("m_", current_month)]

    prediction_grid
}


# %%

## Set target years
years <- 1883:1940

## Final grid that we'll save data to
final_grid <- prediction_grid[c("east", "north")]
final_grid[
    , (ncol(final_grid) + 1):(ncol(final_grid) + nrow(expand.grid(years, 1:12)))
] <- NA

## time_idx is just for tracking what year/month we're on
time_idx <- 1

## Loop through each year
for (current_year in years) {
    year_st <- Sys.time()

    ## Loop through each month
    for (current_month in 1:12) {
        ## Compute the grid
        result <- compute_monthly_grid(
            historic_rain,
            current_year,
            current_month
        )

        ## Add gridded values to output data
        final_grid[, time_idx + 2] <- result$rain
        names(final_grid)[time_idx + 2] <-
            paste0("X", current_year, sprintf("%02d", current_month))

        ## Move on to the next month
        time_idx <- time_idx + 1
    }

    print(paste("year", current_year))
    print(Sys.time() - year_st)
}

# write.csv(
#     final_grid,
#     paste0(
#         "Results/Monthly_Rainfall/Historic_Grids/",
#         "grid_rainfall_", years[1], "_", years[length(years)], ".csv"
#     )
# )
