
################################################################################

# Inverse Distance Weighting - with elevation as a third dimension
# Altitudinal separation between stations can be scaled by a constant value

# From my CV tests this is ineffective, almost always better to
# just do normal IDW

################################################################################

IDW_elevation <- function(y, coords, new_coords,
                          idp = 2, nmax = 8, C_elev = 1) {

    # y                       -observed values of response
    # coords                  -spatial coordinates of observations
    # new_coords              -spatial coordinates to interpolate on to
    #                         coords and new coords should have three columns
    #                         eg (east, north, elev)

    # idp                     -inverse distance weighting power
    # nmax                    -maximum number of neighbours to
    #                          use for prediction
    # C_elev                  -scaling constant for elevation distances

    # Load in packages
    require(sp)
    require(gstat)

    coords <- as.matrix(coords)
    new_coords <- as.matrix(new_coords)

    # For gstat - the data needs to be stored in spatial dataframes
    data <- data.frame(
        y = y, east = coords[, 1],
        north = coords[, 2], elev = coords[, 3]
    )
    newdata <- data.frame(
        east = new_coords[, 1],
        north = new_coords[, 2],
        elev = new_coords[, 3]
    )

    data$elev <- data$elev * C_elev
    newdata$elev <- newdata$elev * C_elev

    coordinates(data) <- c("east", "north", "elev")
    coordinates(newdata) <- c("east", "north", "elev")

    # Interpolate value at each point in newdata
    idw_df <- idw(y ~ 1, data, newdata, idp = idp, nmax = nmax)

    idw_df
}
