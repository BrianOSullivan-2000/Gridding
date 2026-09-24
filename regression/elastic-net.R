
################################################################################

# Regression using elastic-net regularisation
# Elastic-net fits a standard linear regression model, but penalises regression
# parameters that are too high using ridge and lasso penalties

################################################################################

Elastic_Net <- function(df, formula, lambda, alpha = 0.1,
                        seed = 222, nfolds = 10,
                        plot = FALSE, readout = FALSE, ...) {

    # df                      -data frame of response variable and covariates
    # formula                 -regression formula

    # lambda                  -scaling parameter for penalty
    #                         Manually choosing lambda not recommended
    #                         instead find lambda with cv.glmnet in function
    # alpha                   -describes how much of the regression
    #                         is attributed to lasso (alpha = 1) and
    #                         ridge (alpha = 0)
    # seed                    -random seed for cross-validation
    # nfolds                  -number of folds for cross-validation

    # plot                    -plot the output of the regression model
    # readout                 -print a summary of the model

    # Elastic-Net Regularisation is done through the glmnet package
    require(glmnet)
    require(glmnetUtils)

    # Get lambda from CV if unspecified
    if (missing(lambda)) {

        # You should manually set a seed to ensure consistency in CV
        set.seed(seed)
        foldid <- sample(rep(1:nfolds, ceiling(nrow(df) / nfolds)))
        foldid <- head(foldid, nrow(df))

        cvfit <-
            glmnetUtils::cv.glmnet(
                formula, df,
                alpha = alpha,
                use.model.frame = TRUE,
                lambda = 10^seq(log10(5e-7), log10(0.5), length.out = 20),
                foldid = foldid
            )
        lambda <- cvfit$lambda.min
    } else {
        cvfit <- glmnetUtils::glmnet(
            formula, df,
            alpha = alpha, lambda = lambda,
            use.model.frame = TRUE
        )
    }

    # Calculate residuals
    y <- unlist(df[as.character(formula[2])])
    cvfit$residuals <- y - as.vector(predict(cvfit, df, s = lambda))

    # Plot and summary
    if (plot == TRUE) {
        plot(predict(cvfit, df), cvfit$residuals)
        abline(a = 0, b = 0, col = "red")
        qqnorm(cvfit$residuals)
        qqline(cvfit$residuals)
    }
    if (readout == TRUE) {
        print(
            coef(
                glmnet(formula, df,
                    alpha = alpha,
                    lambda = lambda, use.model.frame = TRUE
                )
            )
        )
        print(paste0("Lambda ", lambda))
    }

    cvfit
}