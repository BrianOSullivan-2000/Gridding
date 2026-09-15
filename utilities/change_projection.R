
################################################################################

## Change the projection for spatial coordinates
## This is most commonly used to convert between latitude and longitude
## (EPSG 4326) to the Irish grid (EPSG 29903)

################################################################################

change_projection <- function(coords,
                              source_crs = "EPSG:4326",
                              target_crs = "EPSG:29903") {

    # coords        -spatial coordinates, can be a single coordinate or an array
    #               -it is assumed the coordinates are - (xcoord, ycoord)
    # source_crs    -current projection of coordinates
    # target_crs    -target projection of coordinates

    # Need sf package
    require(sf)

    # Single point
    if (is.numeric(coords) && length(coords) == 2) {
        point <- st_sfc(
            st_point(coords),
            crs = source_crs
        )

        transformed <- st_transform(point, target_crs)
        return(as.numeric(st_coordinates(transformed)))
    }

    # Multiple points (array/matrix)
    coords <- as.matrix(coords)

    if (ncol(coords) != 2) {
        stop("coords must have shape (N, 2)")
    }

    points <- st_as_sf(
        data.frame(x = coords[, 1], y = coords[, 2]),
        coords = c("x", "y"),
        crs = source_crs
    )

    transformed <- st_transform(points, target_crs)

    st_coordinates(transformed)[, 1:2]
}