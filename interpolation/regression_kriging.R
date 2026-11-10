
################################################################################

# Regression Kriging
# Kriging interpolation - the response variable is considered to come from
# a Gaussian Process. Specify a covariance function (eg Matérn) and fit to
# a variogram. Then interpolate using weighted averages.

# For regression kriging, a modelled trend is removed from the data first
# (regression). Then kriging is done on the model residuals. Afterwards, the
# estimated trend is added to the interpolated residuals to get final values

################################################################################

## Load in variogram function
.file <- sys.frame(1)$ofile
.dir <- dirname(normalizePath(.file, mustWork = TRUE))
source(file.path(.dir, "..", "utilities", "variogram.R"))
rm(.file, .dir)

regression_kriging <- function(df, new_df, coords, new_coords,
                               regression_method,
                               nmax = Inf, debug.level = 1, ...) {

    # df                      -data frame of response variable and covariates
    # new_df                  -data frame of covariates for grid locations
    # coords                  -spatial coordinates of observations
    # new_coords              -spatial coordinates to interpolate on to
    #                         coords and new coords should have two columns
    #                         eg (east, north)

    # regression_method       -regression method used for trend (e.g. MLR)

    # nmax                    -maximum number of neighbours to
    #                         use for prediction
    # debug.level             -set to -1 to see progress of kriging
    #                         set to 0 for no printed output

    # Load in packages
    require(sp)
    require(gstat)

    coords <- as.matrix(coords)
    new_coords <- as.matrix(new_coords)

    # Get the residuals
    model <- regression_method(df, coords = coords, ...)
    df$residuals <- model$residuals

    # For gstat - the data needs to be stored in spatial dataframes
    data <- data.frame(
        y = df$residuals,
        east = coords[, 1], north = coords[, 2]
    )
    coordinates(data) <- c("east", "north")

    newdata <- data.frame(
        "east" = new_coords[, 1], "north" = new_coords[, 2]
    )
    coordinates(newdata) <- c("east", "north")

    # Fit spatial variogram to values
    fit_vgm <- spatial_variogram(data, ...)

    # Interpolate values with ordinary kriging
    krig_df <- krige(
        y ~ 1, data, newdata, fit_vgm,
        nmax = nmax,
        debug.level = debug.level
    )
    new_df$residuals <- krig_df$var1.pred

    # Add trend back
    new_df$var1.pred <- predict(model, new_df) + krig_df$var1.pred

    new_df
}