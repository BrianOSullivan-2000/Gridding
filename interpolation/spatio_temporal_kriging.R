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
rm(.file, .dir)

spatio_temporal_kriging <- function(
    y,
    df,
    new_df,
    nmax = Inf,
    debug.level = 1,
    ...
) {
    # y               -observed values of response
    # df              -data frame of input data
    # new_df          -data from of target data
    #                 NOTE both df and new_df require coordinate
    #                 columns east, north, and t

    # nmax            -number of neighbours to use for prediction
    # debug.level     -set to -1 to see progress of kriging
    #                 set to 0 for no printed output

    require(sp)
    require(gstat)
    require(spacetime)

}