
################################################################################

# MLE Kriging
# The model is still just normal ordinary/universal kriging
# the main difference is that the model is estimated using Maximum Likelihood
# Estimation (MLE) rather than minimising least squares from a variogram.

# If the dataset is too large this will be slow,
# but it should be fine for < 2k observations

################################################################################

MLE_kriging <- function(df, new_df, coords, new_coords,
                        init_pars, lower, upper,
                        formula = y ~ 1, cov_function = "Mat", maxit = 1000,
                        useRcpp = FALSE) {

    # df                      -data frame of response variable and covariates
    # new_df                  -data frame of covariates for grid locations
    # coords                  -coordinates of observations (east, north)
    # new_coords              -coordinates to interpolate on to (east, north)
    #                         coords and new coords should have two columns
    #                         eg (east, north)

    # init_pars               -initial covariance parameters
    # lower/upper             -lower and upper limits of covariance parameters
    # formula                 -regression formula (linear regression)
    # cov_function            -type of covariance function to use, can use
    #                         Exponential("Exp"), Matérn("Mat"), Gaussian("Gau")
    # maxit                   -max number of iterations for MLE

    # useRcpp                 -if true, will use faster functions coded in C++

    # Load in necessary packages
    require(sp)
    require(dplyr)
    require(optimParallel)

    ## Setup for optimized code, set up cluster for parallelization,
    ## and source all the Rcpp compiled C++ code on that cluster
    if (useRcpp) {
        solver <- optimParallel
    } else {
        solver <- optim
    }

    # Create spatial objects and get distances between points
    sdf <- df
    new_sdf <- new_df
    sdf$east <- coords[, 1]
    sdf$north <- coords[, 2]
    new_sdf$east <- new_coords[, 1]
    new_sdf$north <- new_coords[, 2]

    coordinates(sdf) <- c("east", "north")
    coordinates(new_sdf) <- c("east", "north")
    ds <- spDists(sdf) / 1000
    cross_ds <- spDists(sdf, new_sdf) / 1000

    # Parameters and covariance function
    if (cov_function %in% c("Exp", "Gau", "Sph")) {
        if (missing(init_pars)) {
            pars = c(1.5, 0.5, 50)
            lower = c(0.001, 0, 1)
            upper = c(Inf, 400, Inf)
        } else if (length(init_pars) != 3) {
            stop("Wrong number of parameters (should be three)")
        }
    }
    if (cov_function == "Mat") {
        if (missing(init_pars)) {
            pars = c(1.5, 0.5, 50, 1.5)
            lower = c(0.001, 0, 1, 0.3)
            upper = c(Inf, Inf, 400, 2.5)
        } else if (length(init_pars) != 4) {
            stop("Wrong number of parameters (should be four)")
        }
    }
    if (!(cov_function %in% c("Exp", "Gau", "Sph", "Mat"))) {
        stop("Model is not available, choose either Exp, Gau, Sph, or Mat")
    }
    if (!missing(init_pars)) {
        pars <- init_pars
    }

    # Response variable and covariates
    y <- model.frame(formula, df)[[1]]
    X <- model.matrix(formula, df)
    new_X <- model.matrix(
        reformulate(
            attr(
                terms(formula),
                "term.labels"
            )
        ), new_df
    )

    # Optimize likelihood numerically to get covariance parameters
    soln <- solver(
        par = pars, fn = spatial_likelihood, ds = ds,
        y = y, X = X, cov_function = cov_function,
        lower = lower, upper = upper,
        control = list(maxit = maxit, parscale = pars, factr = 1e11),
        method = "L-BFGS-B", useRcpp = useRcpp
    )
    pars <- soln$par

    if (useRcpp) {
        new_df$pred <- predict_MVN_cpp(
            pars = pars, ds = ds, cross_ds = cross_ds,
            y = y, X = X, new_X = new_X,
            cov_function = cov_function, beta = numeric(0)
        )
    } else {
        # Get covariance function, inverse, and determinant
        fn <- get(paste0("cov_", cov_function))
        V <- fn(pars, ds)

        # No nugget in cross covariance
        cross_V <- fn(pars, cross_ds)
        Vinv <- solve(V)

        # Solve for regression parameters
        betahat <- solve(t(X) %*% Vinv %*% X) %*% (t(X) %*% (Vinv %*% y))

        new_df$pred <- as.numeric(new_X %*% betahat + t(cross_V) %*%
                                      (Vinv %*% (y - X %*% betahat)))
    }
    new_df
}




spatial_likelihood <- function(pars, ds, y, X, cov_function,
                               beta = NULL, useRcpp = FALSE) {

    if (useRcpp) {
        if (is.null(beta)) {
            beta <- numeric(0)
        }
        negloglik <- spatial_likelihood_cpp(
            pars = pars, ds = ds, y = y, X = X,
            cov_function = cov_function, beta = beta
        )
    } else {
        # Get covariance function, inverse, and determinant
        fn <- get(paste0("cov_", cov_function))
        V <- fn(pars, ds)

        Vchol <- chol(V)
        logdetV <- 2 * sum(log(diag(Vchol)))

        # Get regression parameters
        if (missing(beta)) {
            Vinvy <- backsolve(Vchol, forwardsolve(t(Vchol), y))
            VinvX <- backsolve(Vchol, forwardsolve(t(Vchol), X))

            betahat <- solve(t(X) %*% VinvX) %*% (t(X) %*% Vinvy)
        } else {
            betahat <- beta
        }

        res <- y - (X %*% betahat)
        Vinvres <- backsolve(Vchol, forwardsolve(t(Vchol), res))

        # Calculate negative log-likelihood
        negloglik <- 0.5 * (logdetV + (t(res) %*% Vinvres))
    }

    negloglik
}




## The four different covariance functions available
cov_Exp <- function(pars, ds) {
    sig2 <- pars[1]
    nugget <- pars[2]
    phi <- pars[3]
    V <- sig2 * exp(-ds / phi)
    V[ds == 0] <- V[ds == 0] + nugget
    V
}
cov_Sph <- function(pars, ds) {
    sig2 <- pars[1]
    nugget <- pars[2]
    phi <- pars[3]
    h <- ds / phi
    idx <- h <= 1
    V <- ds * 0

    V[idx] <- sig2 * (1 - 1.5 * h[idx] + 0.5 * h[idx]^3)
    V[ds == 0] <- V[ds == 0] + (sig2 + nugget)
    V
}
cov_Gau <- function(pars, ds) {
    sig2 <- pars[1]
    nugget <- pars[2]
    phi <- pars[3]
    V <- sig2 * exp(-(ds / phi)^2)
    V[ds == 0] <- V[ds == 0] + nugget
    V
}
cov_Mat <- function(pars, ds) {
    sig2 <- pars[1]
    nugget <- pars[2]
    phi <- pars[3]
    kappa <- pars[4]
    h <- ds / phi
    idx <- h > 0
    V <- ds * 0

    V[idx] <- sig2 * ((2^(kappa - 1) * gamma(kappa))^(-1)) *
        ((h[idx])^kappa) *
        besselK(x = h[idx], nu = kappa)
    V[ds == 0] <- V[ds == 0] + (sig2 + nugget)
    V
}