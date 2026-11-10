
################################################################################

# Produce a variogram

# Variograms are useful tools for describing spatial or spatio-temporal
# dependence. A spatial or spatio-temporal covariance function can also be fit
# to them for use in spatial or spatio-temporal interpolation (e.g. kriging)

################################################################################

## Make and fit a spatial variogram

spatial_variogram <- function(
    df,
    flex_vgm = FALSE,
    cutoff = NA, width = NA,

    flex_fit = TRUE,
    vgm_model = "Mat",
    psill = NA, nugget = NA, range = NA,
    kappa = 1,
    fit.method = 7,
    plot_vgm = FALSE,
    ...
) {

    # df               -data (needs response column y)

    # flex_vgm         -automate making empirical variogram (Not recommended)
    # cutoff           -cutoff of empirical variogram
    # width            -width of bins for empirical variogram

    # flex_fit         -automate initial values for theoretical variogram
    #                  (Not recommended)
    # vgm_model        -variogram model
    # psill, nugget,
    # range            -initial partial sill, nugget and range
    # kappa            -kappa value (shape parameter for Matern model)
    # fit.method       -fitting method for vgm
    #                   (see fit.variogram() documentation)
    # plot_vgm         -option to plot the variogram

    # Load in packages
    require(sp)
    require(gstat)

    # Error message if variogram params not
    # provided and flex arguments are False
    if (
        (any(is.na(c(cutoff, width))) && flex_vgm == FALSE) ||
            (any(is.na(c(psill, nugget, range))) && flex_fit == FALSE)
    ) {

        # Find parameters that are missing
        params <- c("cutoff", "width",
                    "psill", "nugget", "range")[is.na(c(cutoff, width,
                                                        psill, nugget, range))]
        params <- paste(params, collapse = ", ")

        warning <- paste0("No values assigned to variogram parameters. 
        Either assign parameters manually (recommended) or use flex_vgm and
        flex_fit arguments. Parameters missing are ", params, ".")

        stop(warning)
    }

    # Flexible way to get cutoff and width of variogram
    if (flex_vgm == TRUE) {
        cutoff <- sqrt((max(df$east) - min(df$east))^2 +
                           (max(df$north) - min(df$north))^2) / 2
        width <- cutoff / 15
    }

    # Empirical variogram
    vgm <- variogram(y ~ 1, data = df, cutoff = cutoff, width = width)

    # Flexible way to get initial variogram parameters
    if (flex_fit == TRUE) {

        nugget <- min(vgm$gamma)

        psill <- median(
            vgm$gamma[
                vgm$dist >= quantile(vgm$dist, 0.75, na.rm = TRUE)
            ],
            na.rm = TRUE
        ) - nugget
        psill <- max(psill, 0)

        range <- vgm$dist[which.min(abs(vgm$gamma - (psill + nugget)))]

    }

    # Set up theoretical variogram - then fit to empirical variogram
    fit_vgm <- vgm(
        model = vgm_model, psill = psill, nugget = nugget,
        range = range, kappa = kappa
    )
    fit_vgm <- fit.variogram(vgm, fit_vgm, fit.method = fit.method)

    if (attr(fit_vgm, "singular")) {

        fit_vgm <- fit.variogram(
            vgm,
            vgm(
                model = vgm_model,
                psill = psill,
                nugget = nugget,
                range = range,
                kappa = kappa
            ),
            fit.sills = c(FALSE, FALSE),
            fit.ranges = TRUE,
            fit.method = fit.method
        )
    }

    if (plot_vgm) {
        print(plot(vgm, fit_vgm))
    }

    fit_vgm
}




spatio_temporal_variogram <- function(
    df,

    flex_vgm = FALSE,
    s_cutoff = NA, s_width = NA,
    t_cutoff = NA, t_width = NA,

    plot_vgm = FALSE
) {

    # df               -data (needs response column y)

    # flex_vgm         -automate making empirical variogram (Not recommended)
    # s_cutoff         -spatial cutoff of empirical variogram
    # s_width          -width of spatial bins for empirical variogram
    # t_cutoff         -temporal cutoff of empirical variogram
    # t_width          -width of temporal bins for empirical variogram

    # plot_vgm         -option to plot the variogram

    fit_vgm <- df
    fit_vgm
}