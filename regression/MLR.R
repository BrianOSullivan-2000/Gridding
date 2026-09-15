
################################################################################

# Multiple linear regression function

################################################################################

MLR <- function(df, formula, step_method = "forward",
                plot = FALSE, readout = FALSE, ...) {

    # df                      -data frame of response variable and covariates
    # formula                 -regression formula
    # step_method             -method to carry out stepwise regression
    #                         (forward, backward, both, none)
    #                         if none, no step-wise regression is calculated

    # plot                    -plot the output of the regression model
    # readout                 -print a summary of the model

    # MASS package required for step-wise
    require(MASS)

    # Standard Regression
    model <- lm(formula = formula, data = df)

    if (tolower(step_method) != "none") {

        # Stepwise Regression
        model <- stepAIC(
            model, method = step_method, trace = FALSE
        )
    }

    # Plot and Summary
    if (plot == TRUE) {
        plot(model)
    }
    if (readout == TRUE) {
        print(summary(model))
    }

    model
}