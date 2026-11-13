
################################################################################

# Produce a variogram

# Variograms are useful tools for describing spatial or spatio-temporal
# dependence. A spatial or spatio-temporal covariance function can also be fit
# to them for use in spatial or spatio-temporal interpolation (e.g. kriging)

################################################################################

## Load in variogram plotting function (for spatio-temporal)
.file <- sys.frame(1)$ofile
.dir <- dirname(normalizePath(.file, mustWork = TRUE))
source(file.path(.dir, "..", "visualisation", "plot_variogram.R"))
rm(.file, .dir)




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

    flex_fit = FALSE,
    s_pars = NA,
    t_pars = NA,
    j_pars = NA,

    vgm_model = c("separable", "sumMetric"),
    spatial_vgm_model = "Exp",
    temporal_vgm_model = "Exp",
    joint_vgm_model = "Exp",
    stAni = NA,
    stAni_range = c(1, 1e8),

    parscales = NA,
    lower_limits = NA,
    upper_limits = NA,
    maxit = 1000,
    fit.method = 7,

    vgm_progress = FALSE,
    plot_vgm = FALSE,
    return_emp_vgm = FALSE,
    no_fit = FALSE
) {

    # df               -input dataframe. Assumes coordinates east, north, t,
    #                  and response variable y

    # flex_vgm         -automate making empirical variogram (Not recommended)
    # s_cutoff         -spatial cutoff of empirical variogram
    # s_width          -width of spatial bins for empirical variogram
    # t_cutoff         -temporal cutoff of empirical variogram
    # t_width          -width of temporal bins for empirical variogram

    # flex_fit         -automatic starting vgm parameters (Not recommended)
    # s_pars           -spatial parameters (sill, range, nugget, kappa)
    # t_pars           -temporal parameters
    # j_pars           -joint parameters (for sum-metric model)

    # vgm_model        -type of spatio-temporal variogram model
    # s_vgm_model      -spatial covariance function
    # t_vgm_model      -temporal covariance function
    # j_vgm_model      -joint covariance function
    # stAni            -spatio-temporal anisotropy constant
    # stAni_range      -search range if fitting spatio-temporal anisotropy

    # parscales        -scales parameters before optimizing over search space
    # lower_limits     -lower limit for fitted parameters
    # upper_limits     -upper limit for fitted parameters
    # maxit            -max number of iterations when fitting variogram
    # fit.method       -fitting method to use for weighting variables
    #                  (see help page for fit.StVariogram)

    # vgm_progress     -show progress bar when generating variogram
    # plot_vgm         -option to plot the variogram
    # return_emp_vgm   -if TRUE, return empirical variogram only
    # no_fit           -if TRUE, don't fit vgm_model just give initial one

    require(sp)
    require(gstat)
    require(spacetime)

    vgm_model <- match.arg(vgm_model)

    if (
        any(is.na(c(s_cutoff, s_width, t_cutoff, t_width))) &&
            flex_vgm == FALSE
    ) {

        ## Warning message
        warning <- paste(
            "No values assigned to empirical variogram parameters.",
            "Either assign them manually (recommended) or use flex_vgm."
        )
        stop(warning)
    }

    if (
        ((any(is.na(c(s_pars, t_pars)))) ||
             (any(is.na(j_pars)) && vgm_model == "sumMetric")) &&
            flex_fit == FALSE
    ) {

        ## Warning message
        warning <- paste(
            "Missing values for variogram model parameters.",
            "Either assign them manually (recommended) or use flex_fit."
        )
        stop(warning)
    }

    ## Heuristics for empirical variogram
    if (flex_vgm) {
        s_cutoff <- sqrt(
            (max(df$east) - min(df$east))^2 + (max(df$north) - min(df$north))^2
        ) / 2
        s_width <- s_cutoff / 15
        t_cutoff <- 12
        t_width <- 1
    }

    ## Create spatial points and dates (with arbitrary start date)
    points <- SpatialPoints(cbind(df$east, df$north))
    dates <- as.Date(df$t, origin = "1970-01-01")

    ## Warning for computation
    if (nrow(df) > 5000) {
        warning(paste(
            "data frame is",
            nrow(df),
            "rows. Consider sampling rows to ease computation"
        ))
    }

    ## Spatio-temporal sparse dataframe
    STIDF <- STIDF(points, dates, data = df)
    STSDF <- as(STIDF, "STSDF")

    ## Empirical variogram
    emp_vgmst <- variogramST(
        y ~ 1,
        data = STSDF,
        tlags = seq(0, t_cutoff, by = t_width),
        cutoff = s_cutoff,
        width = s_width,
        progress = vgm_progress
    )

    if (plot_vgm) {
        zmax <- max(emp_vgmst$gamma, na.rm = TRUE)

        vgm_plot <-
            plot_spatio_temporal_variogram(
                df = df,
                vgmst = emp_vgmst,
                spatial_cutoff = s_cutoff,
                spatial_width = s_width,
                temporal_cutoff = t_cutoff,
                temporal_width = t_width,
                zmax = zmax
            )

        print(vgm_plot)
    }

    if (return_emp_vgm) {
        return(emp_vgmst)
    }

    ## Flexible initial values for variogram
    ## I recommend avoiding this when possible
    if (flex_fit) {

        emp_vgm_spatial <- emp_vgmst[emp_vgmst$timelag == 0, ]
        emp_vgm_temporal <- emp_vgmst[emp_vgmst$spacelag == 0, ]

        ## Only replace a value if it's NA
        replace_parameter_na <- function(parameters, index, fallback) {
            if (length(parameters) < index) {
                length(parameters) <- index
            }
            if (is.na(parameters[index])) {
                parameters[index] <- fallback
            }
            parameters
        }

        ## Nugget
        default_nugget <- min(emp_vgmst$gamma, na.rm = TRUE) / 2
        s_pars <- replace_parameter_na(s_pars, 3, default_nugget)
        t_pars <- replace_parameter_na(t_pars, 3, default_nugget)

        ## Sill
        spatial_sill <- median(
            emp_vgm_spatial$gamma[
                emp_vgm_spatial$dist >= quantile(
                    emp_vgm_spatial$dist, 0.75, na.rm = TRUE
                )
            ],
            na.rm = TRUE
        ) - s_pars[3]
        temporal_sill <- median(
            emp_vgm_temporal$gamma[
                emp_vgm_temporal$timelag >= quantile(
                    emp_vgm_temporal$timelag, 0.75, na.rm = TRUE
                )
            ],
            na.rm = TRUE
        ) - t_pars[3]
        s_pars <- replace_parameter_na(s_pars, 1, spatial_sill)
        t_pars <- replace_parameter_na(t_pars, 1, temporal_sill)

        ## Range
        spatial_range <- emp_vgm_spatial$dist[
            which.min(abs(emp_vgm_spatial$gamma - (s_pars[1] + s_pars[3])))
        ]
        temporal_range <- as.numeric(
            emp_vgm_temporal$timelag[
                which.min(abs(emp_vgm_temporal$gamma - (t_pars[1] + t_pars[3])))
            ]
        )
        s_pars <- replace_parameter_na(s_pars, 2, spatial_range)
        t_pars <- replace_parameter_na(t_pars, 2, temporal_range)

        ## Kappa (for Matérn)
        s_pars <- replace_parameter_na(s_pars, 4, 1)
        t_pars <- replace_parameter_na(t_pars, 4, 1)

        ## For sum-metric model, need joint parameters too
        if (vgm_model == "sumMetric") {
            j_pars <- replace_parameter_na(j_pars, 3, 0)
            j_pars <- replace_parameter_na(
                j_pars,
                1,
                (s_pars[1] + t_pars[1]) / 2
            )
            j_pars <- replace_parameter_na(j_pars, 2, s_pars[2])
            j_pars <- replace_parameter_na(j_pars, 4, 1)
        }
    }

    ## Estimate anisotropy if not provided
    if (is.na(stAni)) {
        stAni <- estiStAni(emp_vgmst, stAni_range)
    }

    if (any(is.na(parscales))) {

        mag <- function(x) {
            10^(ceiling(log10(x))) / 100
        }

        if (vgm_model == "separable") {
            default_parscales <-
                c(
                    mag(s_pars[2]) / 10, mag(s_pars[1]) / 100,
                    mag(t_pars[2]) / 10, mag(t_pars[1]) / 100,
                    mag(s_pars[1] * t_pars[1])
                )

            if (length(parscales) == 5) {
                parscales[is.na(parscales)] <-
                    default_parscales[is.na(parscales)]
            } else {
                parscales <- default_parscales
            }
        }

        if (vgm_model == "sumMetric") {
            default_parscales <-
                c(
                    mag(s_pars[1]), mag(s_pars[2]) / 10, mag(s_pars[3]) / 100,
                    mag(t_pars[1]), mag(t_pars[2]) / 10, mag(t_pars[3]) / 100,
                    mag(j_pars[1]), mag(j_pars[2]) / 10, mag(j_pars[3]) / 100,
                    mag(stAni)
                )

            if (length(parscales) == 10) {
                parscales[is.na(parscales)] <-
                    default_parscales[is.na(parscales)]
            } else {
                parscales <- default_parscales
            }
        }
    }

    ## Heuristics for upper and lower limits
    if (any(is.na(lower_limits)) || any(is.na(upper_limits))) {

        if (vgm_model == "separable") {
            default_lower_limits <-
                c(
                    min(emp_vgmst$dist[emp_vgmst$dist > 0], na.rm = TRUE),
                    0,
                    0.5,
                    0,
                    min(emp_vgmst$gamma, na.rm = TRUE)
                )

            if (length(lower_limits) == 5) {
                lower_limits[is.na(lower_limits)] <-
                    default_lower_limits[is.na(lower_limits)]
            } else {
                lower_limits <- default_lower_limits
            }

            default_upper_limits <-
                c(
                    max(emp_vgmst$dist[emp_vgmst$dist > 0], na.rm = TRUE) * 2,
                    NA,
                    t_pars[2] * 10,
                    NA,
                    max(emp_vgmst$gamma, na.rm = TRUE) * 1.5
                )
            default_upper_limits[c(2, 4)] <-
                default_upper_limits[5] / 3

            if (length(upper_limits) == 5) {
                upper_limits[is.na(upper_limits)] <-
                    default_upper_limits[is.na(upper_limits)]
            } else {
                upper_limits <- default_upper_limits
            }
        }

        if (vgm_model == "sumMetric") {
            default_lower_limits <-
                c(
                    min(emp_vgmst$gamma, na.rm = TRUE) / 2,
                    min(emp_vgmst$dist[emp_vgmst$dist > 0], na.rm = TRUE),
                    0,
                    min(emp_vgmst$gamma, na.rm = TRUE) / 2,
                    0.5,
                    0,
                    min(emp_vgmst$gamma, na.rm = TRUE) / 2,
                    min(emp_vgmst$dist[emp_vgmst$dist > 0], na.rm = TRUE),
                    0,
                    stAni / 100
                )

            if (length(lower_limits) == 10) {
                lower_limits[is.na(lower_limits)] <-
                    default_lower_limits[is.na(lower_limits)]
            } else {
                lower_limits <- default_lower_limits
            }

            emp_vgm_spatial <- emp_vgmst[emp_vgmst$timelag == 0, ]
            emp_vgm_temporal <- emp_vgmst[emp_vgmst$spacelag == 0, ]

            default_upper_limits <-
                c(
                    max(emp_vgm_spatial$gamma, na.rm = TRUE) * 1.5,
                    max(emp_vgmst$dist[emp_vgmst$dist > 0], na.rm = TRUE) * 2,
                    NA,
                    max(emp_vgm_temporal$gamma, na.rm = TRUE) * 1.5,
                    t_pars[2] * 10,
                    NA,
                    max(emp_vgmst$gamma, na.rm = TRUE) * 1.5,
                    min(emp_vgmst$dist[emp_vgmst$dist > 0], na.rm = TRUE),
                    NA,
                    stAni / 100
                )
            default_upper_limits[c(3, 6, 9)] <-
                default_upper_limits[c(1, 4, 7)] / 3

            if (length(upper_limits) == 10) {
                upper_limits[is.na(upper_limits)] <-
                    default_upper_limits[is.na(upper_limits)]
            } else {
                upper_limits <- default_upper_limits
            }
        }
    }

    ## If kappa isn't provided, just put in placeholder
    if (is.na(s_pars[4])) {
        s_pars[4] <- 1
    }
    if (is.na(t_pars[4])) {
        t_pars[4] <- 1
    }

    ## Components of starting variogram
    spatial_vgm <- vgm(
        model = spatial_vgm_model,
        psill = s_pars[1],
        range = s_pars[2],
        nugget = s_pars[3],
        kappa = s_pars[4]
    )
    temporal_vgm <- vgm(
        model = temporal_vgm_model,
        psill = t_pars[1],
        range = t_pars[2],
        nugget = t_pars[3],
        kappa = t_pars[4]
    )

    ## Set up model
    if (vgm_model == "separable") {
        vgm_model_config <- vgmST(
            "separable",
            sill = s_pars[1] * t_pars[1],
            space = spatial_vgm,
            time = temporal_vgm
        )
    }

    if (vgm_model == "sumMetric") {
        if (is.na(j_pars[4])) {
            j_pars[4] <- 1
        }

        ## Components of starting variogram
        joint_vgm <- vgm(
            model = joint_vgm_model,
            psill = j_pars[1],
            range = j_pars[2],
            nugget = j_pars[3],
            kappa = j_pars[4]
        )

        vgm_model_config <- vgmST(
            "sumMetric",
            space = spatial,
            time = temporal,
            joint = joint,
            stAni = stAni
        )
    }

    if (no_fit) {
        fit_vgmst <- vgm_model_config
    } else {

        fit_vgmst <- fit.StVariogram(
            emp_vgmst,
            vgm_model_config,
            fit.method = fit.method,
            stAni = stAni,
            method = "L-BFGS-B",

            control = list(
                parscale = parscales,
                maxit = maxit
            ),
            lower = lower_limits,
            upper = upper_limits
        )
    }

    if (plot_vgm) {
        fitted_surface <- variogramSurface(
            fit_vgmst,
            data.frame(
                spacelag = as.numeric(emp_vgmst$spacelag),
                timelag = as.numeric(emp_vgmst$timelag)
            )
        )

        print(
            plot_spatio_temporal_variogram(
                df = df,
                vgmst = fitted_surface,
                spatial_cutoff = s_cutoff,
                spatial_width = s_width,
                temporal_cutoff = t_cutoff,
                temporal_width = t_width,
                zmax = zmax
            )
        )
    }

    fit_vgmst
}