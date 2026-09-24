
################################################################################

# Random Forest model (RF). Random forest is a regression method that uses an
# ensemble of decision trees when estimating the relationship between the
# response variable and predictors. It's a pretty simple machine learning
# approach, and the version used here is from the ranger package

# Marvin N. Wright and Andreas Ziegler (2017) ranger: A Fast Implementation of
# Random Forests for High Dimensional Data in C++ and R.
# https://doi.org/10.18637/jss.v077.i01

################################################################################

Random_Forest <- function(
    df, formula,
    center_X = TRUE, scale_X = TRUE,

    cv_method = "none", cv_number = 10, cv_index = NULL,
    mtry = 5, splitrule = "variance", min.node.size = 5,
    seed = 222, plot = FALSE, readout = FALSE,
    ...
) {

    # df                      -data frame of response variable and covariates
    # formula                 -regression formula
    # center_X                -if TRUE, centre X (recommended)
    # scale_X                 -if TRUE, scale the predictors (recommended)

    # cv_method               -CV method used (if any) to tune parameters
    # cv_number               -describes number of folds/samples used for CV
    # cv_index                -you can manually set indices for folds

    # mtry                    -number of predictors at each split in the RF
    # splitrule               -by what metric splitting decisions are made
    # min.node.size           -min number of observations allowed from a node

    # seed                    -random seed for reproducibility
    # plot                    -plot the output of the regression model
    #                         for RF, this is the significance of each predictor
    # readout                 -print a summary of the model

    require(caret)
    require(ranger)

    ## Get design matrix and response data in vector form
    X <- as.matrix(model.matrix(formula, df)[, -1])
    y <- as.numeric(model.response(model.frame(formula, df)))

    ## Tuning parameters when training random forest model
    tune_grid <- expand.grid(
        mtry = mtry,
        splitrule = splitrule,
        min.node.size = min.node.size
    )

    ## Train the model
    ## Note - train() from the caret package let's you do a much more involved
    ## model training than what is done here. I want to keep it simple though

    set.seed(seed)
    model <- train(
        X, y,
        method = "ranger",
        trControl = trainControl(
            method = cv_method,
            number = cv_number,
            index = cv_index
        ),
        tuneGrid = tune_grid,
        preProcess = c("center", "scale")[c(center_X, scale_X)],
        num.trees = 500,
        importance = ifelse(cv_method == "none", "none", "permutation"),
        verbose = TRUE
    )

    ## Get model residuals
    model$residuals <- y - predict(model, newdata = X)

    # Plot and Summary
    if (plot == TRUE) {
        plot(model)
    }
    if (readout == TRUE) {
        print(summary(model))
    }

    model
}