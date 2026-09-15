
################################################################################

# Plot a nice grid over the island

################################################################################

base_grid_plot <- function(grid, island_outline, stations = NA,
                           response = "y", xy_names = c("easting", "northing"),

                           title = "", legend_title = "", yunit = NA,

                           grid_size = 0.1, station_size = 2,

                           colour_palette = c("#ffffd9", "#1d91c0", "#081d58"),
                           breaks = NA, labels = NA,
                           nbreaks = 20, bias = 1,
                           ymin = min(grid[response]),
                           ymax = max(grid[response]),
                           plot_destination = NA) {

    # grid                  -dataframe with response value and coordinates
    #                       coordinates need to be labelled east, north
    # island_outline        -outline of the island for the plot
    #                       (x,y coordinates)
    # stations              -stations used to make grid
    # response              -string identifying variable to plot
    # xy_names              -label for plotting coordinates
    #                       eg east/north, lon/lat

    # title                 -title of plot
    # legend_title          -title of legend
    # yunit                 -unit of plot variable

    # grid_size             -size of grid points
    # station_size          -size of stations (if any)

    # colour_palette        -palette to use for plot
    # breaks, labels        -can add breaks and labels yourself
    # nbreaks               -number of bins for colour palette
    # bias                  -bias argument for colorRampPalette (see docs)
    # ymin, ymax            -limits of colour palette (I recommend doing
    #                       this manually for nice round numbers)
    # plot_destination      -can save plot

    # Load in packages
    require(ggplot2)
    require(scales)

    # Have a consistent name for the response variable & xy coordinates
    grid$yplot <- unlist(grid[response])
    names(grid)[names(grid) == xy_names[1]] <- "east"
    names(grid)[names(grid) == xy_names[2]] <- "north"

    if (any(!is.na(stations))) {
        stations$yplot <- unlist(stations[response])
        names(stations)[names(stations) == xy_names[1]] <- "east"
        names(stations)[names(stations) == xy_names[2]] <- "north"
    }

    # Set up the colour palette for the plot
    if ((nbreaks > 0) && any(is.na(breaks))) {
        breaks <- pretty(c(ymin, ymax), n = nbreaks, min.n = nbreaks)
        nbreaks <- length(breaks)
        colour_palette <-
            colorRampPalette(colour_palette, bias = bias)(nbreaks)
    } else if (any(is.na(breaks))) {
        colour_palette <-
            colorRampPalette(colour_palette, bias = bias)
    } else {
        nbreaks <- length(breaks)
        colour_palette <-
            colorRampPalette(colour_palette, bias = bias)(nbreaks)
    }

    # Palettes for binned vs. continuous
    if (!any(is.na(breaks))) {

        if (any(is.na(labels))) {
            # Add a nice unit
            if (!is.na(yunit)) {
                labels <- rep("", nbreaks)

                idx <- round(seq(1, nbreaks,
                    length.out = min(10, nbreaks)
                ))
                labels[idx] <- paste(as.character(breaks[idx]), yunit)
            } else {
                labels <- rep("", nbreaks)

                idx <- round(seq(1, nbreaks,
                    length.out = min(10, nbreaks)
                ))
                labels[idx] <- as.character(breaks[idx])
            }
        }

        # color for grid - fill for stations
        color_scale <-
            scale_colour_stepsn(
                colours = colour_palette,
                values = scales::rescale(breaks, to = c(0, 1)),
                breaks = breaks, labels = labels,
                limits = c(min(breaks), max(breaks))
            )

        if (any(!is.na(stations))) {
            fill_scale <-
                scale_fill_stepsn(
                    colours = colour_palette,
                    values = scales::rescale(breaks, to = c(0, 1)),
                    breaks = breaks, labels = labels,
                    limits = c(min(breaks), max(breaks))
                )
        }
    } else {

        # color for grid - fill for stations
        color_scale <-
            scale_colour_gradientn(
                colours = colour_palette(256),
                limits = c(signif(ymin, 3), signif(ymax, 3))
            )

        if (any(!is.na(stations))) {
            fill_scale <-
                scale_fill_gradientn(
                    colours = colour_palette(256),
                    limits = c(signif(ymin, 3), signif(ymax, 3))
                )
        }
    }

    # Plotting
    ireland_plot <-
        ggplot(grid) +

        # Title
        ggtitle(title) +

        # Border
        theme_void() +
        geom_point(
            data = island_outline, aes(x, y),
            colour = "gray50", size = 0.5
        ) +

        # Aesthetics
        theme(
            aspect.ratio = 1.3,
            text = element_text(size = 11),
            plot.title = element_text(margin = margin(b = 2)),
            axis.title = element_blank(),
            panel.grid.major = element_line(
                color = gray(.5), linetype = "dotted", linewidth = 0.2
            ),
            panel.grid.minor = element_line(
                color = gray(.5), linetype = "dotted", linewidth = 0.2
            ),
            panel.background = element_rect(fill = "transparent")
        ) +
        labs(colour = legend_title) +
        guides(fill = "none") +

        theme(legend.key.size = unit(1, "cm"),
              legend.key.width = unit(0.5, "cm"),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              plot.title = element_text(size = 9.5, hjust = 0.7,
                                        vjust = 0, face = "bold"),
              plot.subtitle = element_text(size = 9, hjust = 0.5),
              plot.margin = margin(0, 3, 0, 0, "mm")) +

        # Gridded values
        geom_point(data = grid, aes(x = east, y = north, color = yplot),
                   size = grid_size, alpha = 1, show.legend = TRUE) +

        color_scale

    # Plot the stations too if they're included
    if (any(!is.na(stations))) {
        # Stations
        ireland_plot <-
            ireland_plot +
            geom_point(
                data = stations, aes(x = east, y = north, fill = yplot),
                size = station_size, shape = 21, show.legend = TRUE
            ) +
            fill_scale
    }

    if (!is.na(plot_destination)) {
        ggsave(plot_destination, ireland_plot)
    }

    ireland_plot
}