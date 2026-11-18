################################################################################

# Ordinary Kriging (spatio-temporal)
# Kriging interpolation - the response variable is considered to come from
# a Gaussian Process. Specify a covariance function and fit to a variogram,
# then interpolate using weighted averages.

# Moving from a spatial to a spatio-temporal context for kriging isn't trivial
# You have to specify a spatio-temporal basis function and produce a spatio-
# temporal variogram. From there, you do kriging as normal

################################################################################

## Load in variogram function
.file <- sys.frame(1)$ofile
.dir <- dirname(normalizePath(.file, mustWork = TRUE))
source(file.path(.dir, "..", "utilities", "variogram.R"))

spatio_temporal_kriging <- function(
    df,
    new_df,
    response = "y",
    nmax = Inf,
    progress = FALSE,
    ...
) {
    # df              -data frame of input data
    # new_df          -data from of target data
    #                 NOTE both df and new_df require coordinate
    #                 columns east, north, and t
    # response        -target variable to be interpolated

    # nmax            -number of neighbours to use for prediction
    # progress        -if TRUE, show progress bar during kriging

    require(sp)
    require(gstat)
    require(spacetime)

    ## Get response data
    y <- df[[response]]
    df$y <- y

    ## Get fitted variogram
    vgmst <- spatio_temporal_variogram(
        df, ...
    )
    stAni <- vgmst$stAni
    vgmst <- vgmst$vgmst

    ## Get spatial points
    points <- SpatialPoints(
        cbind(df$east, df$north)
    )
    ## Create spatio-temporal sparse dataframes
    STIDF <- STIDF(
        points,
        as.Date(df$t, origin = "1970-01-01"),
        data = df
    )
    ST <- as(STIDF, "STSDF")

    ## Repeat for new data
    new_points <- SpatialPoints(
        cbind(new_df$east, new_df$north)
    )
    new_STIDF <- STIDF(
        new_points,
        as.Date(new_df$t, origin = "1970-01-01"),
        data = new_df
    )
    new_ST <- as(new_STIDF, "STSDF")

    ## Kriging
    krig_df <- krigeST(
        y ~ 1,
        ST,
        new_ST,
        vgmst,
        stAni = stAni,
        nmax = nmax,
        progress = progress
    )

    krig_df
}