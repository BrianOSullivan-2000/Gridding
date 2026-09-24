
################################################################################

# Extreme Gradient Boosting model (XGBoost). XGBoost is a regression method that
# uses an ensemble of decision trees, where trees are built sequentially and
# each new tree is trained to improve the predictions of the previous trees.
# It's a flexible machine learning approach, and the version used here is from
# the xgboost package

# Tianqi Chen and Carlos Guestrin (2016) XGBoost: A Scalable Tree Boosting
# System.
# https://doi.org/10.1145/2939672.2939785

################################################################################

XGBOOST <- function(
    df, formula,
    center_X = FALSE, scale_X = FALSE,

    coords = NULL,

    nrounds = 500,
    max_depth = 6,
    learning_rate = 0.1,
    min_child_weight = 1,
    subsample = 1,
    colsample_bytree = 1,

    seed = 222, plot = FALSE, readout = FALSE,
    ...
) {

    # df                  -data frame of response variable and covariates
    # formula             -regression formula
    # center_X            -if TRUE, centre X
    # scale_X             -if TRUE, scale the predictors

    # coords               -included to avoid naming conflict

    # nrounds              -number of boosting rounds (trees/epochs)
    # max_depth            -maximum depth of each tree (overall complexity)
    # learning_rate        -step size shrinkage used in boosting (eta)
    # min_child_weight     -minimum sum of observation weights required
    # subsample            -proportion of observations sampled for each tree
    # colsample_bytree     -proportion of predictors sampled for each tree
    #                      Alongside subsample, very useful for regularization

    # seed                -random seed for reproducibility
    # plot                -plot the output of the model
    # readout             -print a summary of the model

    require(xgboost)

    ## Get design matrix and response
    X <- model.matrix(
        delete.response(terms(formula)),
        data = df
    )
    ## Remove intercept
    X <- X[, -1, drop = FALSE]

    X <- matrix(
        as.numeric(X),
        nrow = nrow(X),
        ncol = ncol(X),
        dimnames = dimnames(X)
    )

    y <- as.numeric(
        model.response(
            model.frame(formula, df)
        )
    )

    X <- scale(
        X,
        center = center_X,
        scale = scale_X
    )
    X_center <- attr(X, "scaled:center")
    X_scale <- attr(X, "scaled:scale")

    ## Train XGBoost model
    set.seed(seed)

    xgb_model <- xgboost(
        x = X,
        y = y,
        objective = "reg:squarederror",
        nrounds = nrounds,
        max_depth = max_depth,
        learning_rate = learning_rate,
        min_child_weight = min_child_weight,
        subsample = subsample,
        colsample_bytree = colsample_bytree,
        ...
    )

    ## Get model residuals
    residuals <- y - predict(xgb_model, newdata = X)

    ## Create wrapper model (for prediction)
    model <- list(
        model = xgb_model,
        residuals = residuals,
        formula = formula,
        center_X = center_X,
        scale_X = scale_X,
        X_center = X_center,
        X_scale = X_scale
    )

    class(model) <- "XGBOOST"

    ## Plot
    if (plot == TRUE) {
        xgb.plot.importance(
            xgb.importance(model = xgb_model)
        )
    }

    ## Readout
    if (readout == TRUE) {
        print(xgb_model)
    }

    model
}




predict.XGBOOST <- function(object, newdata, ...) {

    X_new <- model.matrix(
        delete.response(terms(object$formula)),
        data = newdata
    )
    X_new <- X_new[, -1, drop = FALSE]
    X_new <- matrix(
        as.numeric(X_new),
        nrow = nrow(X_new),
        ncol = ncol(X_new),
        dimnames = dimnames(X_new)
    )

    if (object$center_X) {
        X_new <- sweep(
            X_new,
            2,
            object$X_center,
            FUN = "-"
        )
    }

    if (object$scale_X) {
        X_new <- sweep(
            X_new,
            2,
            object$X_scale,
            FUN = "/"
        )
    }

    predict(
        object$model,
        newdata = X_new,
        ...
    )
}