
################################################################################

# Generalised Least Squares Regression
# Calculate a covariance between locations using a variogram, then use that
# when fitting regression parameters

################################################################################

GLS <- function(
    df, formula, fit_vgm,
    multi_GLS = FALSE, n_iter = 5,
    plot = FALSE, readout = FALSE, ...
) {

    # df                      -data frame of response variable and covariates
    # formula                 -regression formula

    # fit_vgm                 -variogram (covariance structure of process)
    # multi_GLS               -recursively fit the variogram several times
    # n_iter                  -number of recursions of multi_GLS
    #                         you don't need many (5 is definitely enough)

    # plot                    -plot the output of the regression model
    # readout                 -print a summary of the model

    # Load in necessary packages
    require(sp)
    require(gstat)

    # Create spatial object and get distances between points
    sdf <- df
    coordinates(sdf) <- c("east", "north")
    ds <- spDists(sdf)

    # Use this if a variogram is provided
    if (!missing(fit_vgm)) {
        # Using fitted variogram - get covariance between all observations
        V <- variogramLine(fit_vgm, dist_vector = ds, covariance = TRUE)

        # Actual formula for GLS
        X <- model.matrix(formula, df)
        y <- unlist(df[as.character(formula[2])])
        beta <- solve((t(X) %*% solve(V, X)), (t(X) %*% (solve(V, y))))

        # Set up a model object (just a way to store all relevant info)
        model <- list()
        class(model) <- "gls"
        model$formula <- formula
        model$coefficients <- beta
        model$fitted.values <- as.numeric(X %*% beta)
        model$residuals <- y - model$fitted.values
    } else {

        # Make and fit a variogram if you don't already have one
        model <- lm(formula = formula, data = df)
        sdf$residuals <- model$residuals
        fit_vgm <- spatial_variogram_GLS(sdf, ...)

        iter <- ifelse(multi_GLS, n_iter, 1)
        for (i in 1:iter) {

            if (i > 1) {
                sdf$residuals <- model$residuals
                fit_vgm <- spatial_variogram_GLS(sdf, ...)
            }

            # Using fitted variogram - get covariance between all observations
            V <- variogramLine(fit_vgm, dist_vector = ds, covariance = TRUE)

            # Actual formula for GLS
            X <- model.matrix(formula, df)
            y <- unlist(df[as.character(formula[2])])
            beta <- solve((t(X) %*% solve(V, X)), (t(X) %*% (solve(V, y))))

            # Set up a model object (just a way to store all relevant info)
            model <- list()
            class(model) <- "gls"
            model$formula <- formula
            model$coefficients <- beta
            model$fitted.values <- as.numeric(X %*% beta)
            model$residuals <- y - model$fitted.values

        }
    }

    # Plot and summary
    if (plot == TRUE) {
        plot(
            model$fitted.values, model$residuals, ,
            xlab = "Fitted values", ylab = "Residuals"
        )
        abline(a = 0, b = 0, col = "red")
        qqnorm(model$residuals)
        qqline(model$residuals)
    }
    if (readout == TRUE) {
        print(paste0("Fitted coefficients for ", deparse(formula), ":"))
        print(signif(model$coefficients, 3))
    }

    model
}




## Make and fit a spatial variogram

spatial_variogram_GLS <- function(df, flex_vgm = FALSE, flex_fit = FALSE,
                                  cutoff = NA, width = NA,
                                  vgm_model = "Mat",
                                  psill = NA, nugget = NA, range = NA,
                                  kappa = 1, plot_vgm = FALSE, ...) {

    # df               -data
    # flex_vgm         -automate making empirical variogram (Not recommended)
    # flex_fit         -automate initial values for theoretical variogram
    #                  (Not recommended)

    # cutoff           -cutoff of empirical variogram
    # width            -width of bins for empirical variogram

    # vgm_model        -variogram model
    # psill, nugget,
    # range            -initial partial sill, nugget and range
    # kappa            -kappa value (shape parameter for Matern model)
    # plot_vgm         -option to plot the variogram

    # Load in packages
    require(sp)
    require(gstat)

    # Error message if variogram params not provided and
    # flex arguments are False
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
    f <- as.formula("residuals ~ 1")
    vgm <- variogram(f, data = df, cutoff = cutoff, width = width)

    # Flexible way to get initial variogram parameters
    if (flex_fit == TRUE) {
        nugget <- mean(vgm$gamma[1:3])
        psill <- mean(vgm$gamma) - nugget
        range <- vgm$dist[which.min(abs(vgm$gamma - (psill + nugget)))]
    }

    # Set up theoretical variogram - then fit to empirical variogram
    fit_vgm <- vgm(
        model = vgm_model, psill = psill, nugget = nugget,
        range = range, kappa = kappa
    )
    fit_vgm <- fit.variogram(vgm, fit_vgm)

    if (plot_vgm) {
        print(plot(vgm, fit_vgm))
    }

    fit_vgm
}




predict.gls <- function(gls, df) {
    X <- model.matrix(delete.response(terms(gls$formula)), df)
    y <- as.numeric(X %*% gls$coefficients)
    y
}
