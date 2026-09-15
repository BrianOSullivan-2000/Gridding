
################################################################################

# Principal Component Regression
# Take the covariates (east, north, elev) and combine the ones that have high
# covariance into principal components. These components are used for regression

################################################################################


PCR <- function(df, formula, center = TRUE, scale = TRUE,
                plot = FALSE, readout = FALSE, elbow_plot = FALSE,
                validation = "CV",
                ncomp = NA, plsr = FALSE, ...) {

    # df              -dataframe (response variable and explanatory variables)
    # formula         -regression formula
    # center          -center the data around a mean before making components
    # scale           -scale the data (standardise)
    #                 highly recommended to set to True

    # plot            -plot the output of the regression model
    # readout         -print a summary of the model
    # elbow_plot      -create an elbow plot of the data
    # validation      -include validation method for selecting #components
    # plsr            -can opt to use plsr instead of pcr
    #                  (partial least squares)

    # Done with pls package
    require(pls)

    # Create either PLSR or PCR regression model
    if (plsr == TRUE) {
        model <- plsr(
            formula = formula, data = df,
            center = center, scale = scale, validation = validation
        )
    } else {
        model <- pcr(
            formula = formula, data = df,
            center = center, scale = scale, validation = validation
        )
    }

    # Find optimal number of components - option to output an elbow plot
    ncomp <- selectNcomp(model, plot = elbow_plot)

    if (plsr == TRUE) {
        model <- plsr(
            formula = formula, data = df,
            center = center, scale = scale,
            validation = validation, ncomp = ncomp
        )
    } else {
        model <- pcr(
            formula = formula, data = df,
            center = center, scale = scale,
            validation = validation, ncomp = ncomp
        )
    }

    # Plotting and summary
    if (plot == TRUE) {
        plot(predict(model, df), model$residuals)
        abline(a = 0, b = 0, col = "red")
        qqnorm(model$residuals)
        qqline(model$residuals)
    }

    if (readout == TRUE) {
        print(model$coefficients)
    }

    model$residuals <- model$residuals[, , ncomp]
    model$formula <- formula
    class(model) <- "pcr"

    model
}




# I've had to manually recode the predict function
# in the pls package they predict for every ncomp

predict.pcr <- function(pcr, df) {

    ncomp <- pcr$ncomp

    # Get model matrix
    Terms <- delete.response(terms(pcr))
    X <- model.matrix(Terms, df)
    X <- X[, !(colnames(X) %in% "(Intercept)"), drop = FALSE]
    X <- scale(X, center = pcr$Xmeans, scale = pcr$scale)

    # Get coefficients and intercept
    B <- coef(pcr, ncomp = ncomp)[, , ncomp]
    intercept <- as.numeric(pcr$Ymean - crossprod(pcr$Xmeans, B))

    # Predict
    y_pred <- X %*% B + intercept
    y_pred
}
