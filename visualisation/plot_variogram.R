
################################################################################

# Variograms are great tools for exploratory data analysis
# This script can plot spatial and spatio-temporal variograms to help
# investigate spatial and spatio-temporal dependence in data, respectively

################################################################################

plot_spatial_variogram <- function(
    df,
    response = "y",
    cutoff = 400000, width = 20000,
    title = "",
    plot_curve = FALSE, fit_trend = FALSE,
    psill = NA, nugget = NA, range = NA, kappa = NA
) {
    # df               -data (needs coordinate columns east + north)
    # response         -name of quantity to be plotted
    # cutoff           -cutoff of empirical variogram
    # width            -width of bins for empirical variogram
    # title            -title of plot
    # plot_curve       -if TRUE, fit a curve (Matérn function) to the data
    # fit_trend        -if TRUE, will fit curve to variogram
    # psill, nugget,
    # range, kappa     -initial partial sill, nugget, range, and kappa

    ## Load in packages
    require(sp)
    require(gstat)
    require(ggplot2)
    require(ggtext)

    ## Empirical variogram
    f <- as.formula(paste0(response, " ~ 1"))
    coordinates(df) <- c("east", "north")
    vgm <- variogram(f, data = df, cutoff = cutoff, width = width)

    ## Flexible way to get initial variogram parameters
    if (plot_curve == TRUE) {
        if (any(is.na(c(psill, nugget, range, kappa)))) {
            message("No covariance params provided, using flexible params")
            fit_trend <- TRUE

            ## Set nugget as min value
            nugget <- min(vgm$gamma)

            ## Sill defined using 75th quantile
            psill <- median(
                vgm$gamma[
                    vgm$dist >= quantile(vgm$dist, 0.75, na.rm = TRUE)
                ],
                na.rm = TRUE
            ) - nugget
            psill <- max(psill, 0)

            ## Range define as where sill is met
            range <- vgm$dist[which.min(abs(vgm$gamma - (psill + nugget)))]

            ## Default to Exponential
            kappa <- 0.5
        }

        ## Define Matérn function
        fit_vgm <- vgm(
            model = "Mat", psill = psill, nugget = nugget,
            range = range, kappa = kappa
        )

        ## Fit to empirical variogram
        if (fit_trend) {
            fit_vgm <- fit.variogram(vgm, fit_vgm)
        }
    }

    ## Plotting
    vgm_plot <- ggplot(vgm) +
        geom_point(
            mapping = aes(x = dist, y = gamma)
        ) +
        xlab("Distance (m)") +
        ylab("Semivariance") +
        ylim(c(0, max(vgm$gamma, na.rm = TRUE))) +
        scale_x_continuous(
            labels = scales::label_number()
        ) +
        ggtitle(title)

    ## Overlay curve
    if (plot_curve) {

        curve_vals <- variogramLine(
            fit_vgm,
            dist_vector = seq(1e-6, cutoff, length.out = 2000)
        )
        vgm_plot$layers <- append(
            vgm_plot$layers,
            list(geom_line(
                data = curve_vals,
                aes(x = dist, y = gamma),
                color = "red"
            )),
            after = 0
        )

        ## Annotate covariance parameters
        label_df <- data.frame(
            x = Inf,
            y = -Inf,
            label = paste0(
                "&sigma;<sup>2</sup> = ", fit_vgm[2, 2], "<br>",
                "&tau;<sup>2</sup> = ", fit_vgm[1, 2], "<br>",
                "&psi; = ", fit_vgm[2, 3], "<br>",
                "&kappa; = ", fit_vgm[2, 4]
            )
        )
        vgm_plot +
            geom_textbox(
                data = label_df,
                aes(x = x, y = y, label = label),
                hjust = 1.05,
                vjust = -0.05,
                halign = 0,
                width = grid::unit(1, "in"),
                lineheight = 1.5,
                fill = "white",
                box.colour = "black"
            )
    }

    vgm_plot
}




plot_spatio_temporal_variogram <- function(
    df,
    response = "y",
    spatial_cutoff = 400000, spatial_width = 20000,
    temporal_cutoff = 12, temporal_width = 1,
    time_label = "Time"
) {

    ## TODO - Add in some curve specification/fitting
    ## it's much harder to get heuristic starting values
    ## for this compared to the spatial-only case

    # df               -data (needs coordinate columns east + north + t)
    # response         -name of quantity to be plotted
    # spatial_cutoff   -spatial cutoff of empirical variogram
    # spatial_width    -spatial width of empirical variogram
    # temporal_cutoff  -temporal cutoff of empirical variogram
    # temporal_width   -temporal width of empirical variogram
    # time_label       -label for temporal axis (e.g. "Days")

    ## Load in packages
    require(sp)
    require(gstat)
    require(lattice)

    ## Empirical variogram formula
    f <- as.formula(paste0(response, " ~ 1"))

    ## Need to set a start date (arbitrary)
    start_date <- "1970-01-01"

    ## Create spatio-temporal sparse dataframe
    points <- SpatialPoints(cbind(df$east, df$north))
    times <- as.Date(df$t, origin = start_date)
    STIDF <- STIDF(points, times, data = df)
    STSDF <- as(STIDF, "STSDF")

    ## Calculate variogram
    vgmst <- variogramST(
        f,
        data = STSDF,
        tlags = seq(0, temporal_cutoff, by = temporal_width),
        cutoff = spatial_cutoff, width = spatial_width,
        progress = TRUE
    )

    ## Need to make space and time lags numeric for plotting
    spacelags <- seq(0, spatial_cutoff, by = spatial_width)
    timelags <- seq(0, temporal_cutoff, by = temporal_width)
    vgmst$spacelag <- rep(spacelags, length(timelags))
    vgmst$timelag <- rep(timelags, each = length(spacelags))
    vgmst$spacelag <- as.numeric(vgmst$spacelag)
    vgmst$timelag <- as.numeric(vgmst$timelag)

    ## Plotting
    vgmst_plot <- wireframe(
        gamma ~ spacelag * timelag,
        vgmst,
        drape = TRUE,
        col.regions = col,
        colorkey = FALSE,
        as.table = TRUE,
        scales = list(arrows = FALSE, cex = 1),
        par.settings = list(
            axis.text = list(cex = 1),
            axis.line = list(col = "transparent")
        ),
        xlab = list("Distance (m)", cex = 1.5, rot = 30),
        ylab = list(time_label, cex = 1.5, rot = -35),
        zlab = "", zlim = c(0, max(vgmst$gamma)), zoom = 0.95
    )

    vgmst_plot
}
