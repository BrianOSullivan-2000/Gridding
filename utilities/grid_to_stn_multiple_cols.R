
#-----------------------------------------------------#
## AUTHORED BY BARRY COONAN - 18th November 2026 ##
#-----------------------------------------------------#

## Calculate weighted values for a location from a grid

grid_to_stn_multiple_cols <- function(station_data, grid) {

    en_names <- 1 # identify the coordinate label types

    # station coordinates
    if (colnames(station_data)[2] == "easting") {
        en_names <-  2
        colnames(station_data)[2:3] <- c("east", "north")
    }

    #  station_data <- c(station_data$east, station_data$north)
    #  colnames(grid)[3] <- "param"

    for (stn in seq_along(nrow(station_data))) {
        # if point of interest not on a grid line or is not a grid point
        # then use 4 nearest points
        near_points <- matrix(
            c(
                round(floor(station_data$east[stn] / 1000), 0) * 1000,
                round(ceiling(station_data$east[stn] / 1000), 0) * 1000,
                round(floor(station_data$east[stn] / 1000), 0) * 1000,
                round(ceiling(station_data$east[stn] / 1000), 0) * 1000,
                round(floor(station_data$north[stn] / 1000), 0) * 1000,
                round(floor(station_data$north[stn] / 1000), 0) * 1000,
                round(ceiling(station_data$north[stn] / 1000), 0) * 1000,
                round(ceiling(station_data$north[stn] / 1000), 0) * 1000
            ),
            ncol = 2, nrow = 4
        )

        # if location on a grid point or between 2 grid points
        near_points <- near_points[!duplicated(near_points), ]

        near_points <- as.data.frame(near_points)

        # get nearest grid point weights
        # (only if the location is not on a grid point)
        if (nrow(near_points) > 1) {
            # weights => get distance from target point to grid points
            # weight = inverse of distance squared

            # point of interest not a grid point or on gridline
            if (nrow(near_points) == 4) {
                point_weights <- c(
                    1 / ((station_data$east[stn] - near_points[1, 1])^2 +
                             (station_data$north[stn] - near_points[1, 2])^2),
                    1 / ((station_data$east[stn] - near_points[2, 1])^2 +
                             (station_data$north[stn] - near_points[2, 2])^2),
                    1 / ((station_data$east[stn] - near_points[3, 1])^2 +
                             (station_data$north[stn] - near_points[3, 2])^2),
                    1 / ((station_data$east[stn] - near_points[4, 1])^2 +
                             (station_data$north[stn] - near_points[4, 2])^2)
                )
                near_points$weight <- point_weights
                near_points$weight <-
                    near_points$weight / sum(near_points$weight)

            }

            # point of interest on a gridline
            if (nrow(near_points) == 2) {
                point_weights <- c(
                    ((station_data$east[stn] - near_points[1, 1])^2 +
                         (station_data$north[stn] - near_points[1, 2])^2),
                    ((station_data$east[stn] - near_points[2, 1])^2 +
                         (station_data$north[stn] - near_points[2, 2])^2)
                )
                near_points$weight <- point_weights
                near_points$weight <-
                    near_points$weight / sum(near_points$weight)

            }
        }

        colnames(near_points)[1:2] <- c("east", "north")

        # combine the geo-data and weights data
        near_points_data <- merge(
            near_points, grid,
            by = c("east", "north"), all.x = TRUE
        )
        near_points_data <- na.omit(near_points_data)

        if (nrow(near_points_data) < 4) {
            near_points_data$weight <-
                near_points_data$weight / sum(near_points_data$weight)
        }

        weighted_params <- apply(
            near_points_data$weight *
                near_points_data[, 4:ncol(near_points_data)],
            2, sum
        )

        # combine with date and station info
        near_points_data <- cbind(
            station_data$stno[stn],
            station_data$east[stn],
            station_data$north[stn],
            t(apply(
                near_points_data$weight *
                    qnear_points_data[, 4:ncol(near_points_data)],
                2, sum
            ))
        )

        if (stn == 1) {
            stn_values <- near_points_data
        } else {
            stn_values <- rbind(stn_values, near_points_data)
        }
    }

    stn_values <- data.frame(stn_values)

    colnames(stn_values)[1:3] <- colnames(station_data)

    stn_values$stno <- as.integer(stn_values$stno)
    stn_values$east <- as.integer(stn_values$east)
    stn_values$north <- as.integer(stn_values$north)
    stn_values[, 4:ncol(stn_values)] <-
        round(stn_values[, 4:ncol(stn_values)], 1)

    if (en_names == 2) colnames(stn_values)[2:3] <- c("easting", "northing")

    stn_values

}
