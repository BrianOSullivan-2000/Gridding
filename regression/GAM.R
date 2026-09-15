
################################################################################

# Generalized Additive Model (GAM)

################################################################################

GAM <- function(df, formula,
                # regression_estimation = "likelihood",

                plot = FALSE, readout = FALSE, ...) {

    # df                      -data frame of response variable and covariates
    # formula                 -regression formula
    # regression_estimation   -fitting method (likelihood or Bayesian)

    # plot                    -plot the output of the regression model
    # readout                 -print a summary of the model

    require(mgcv)
    model <- gam(formula, data = df, ...)

    # Likelihood estimation
    # if(regression_estimation == "likelihood"){
    #   require(mgcv)
    #   model <- gam(formula, data = df, ...)
    # }

    # # Bayesian estimation
    # else if(regression_estimation == "Bayesian"){
    #   require(bamlss)
    #   model <- bamlss(formula, data = df, ...)
    #   model$residuals <- residuals(model, type = "response")
    #
    #   if(model$family$family == "cnorm"){
    #     class(model) <- c("GAM_bamlss_cnorm", class(model))
    #   }
    #   else{
    #     class(model) <- c("GAM_bamlss", class(model))
    #   }
    # }

    # Plot and Summary
    if (plot == TRUE) {
        plot(model)
    }
    if (readout == TRUE) {
        print(summary(model))
    }

    model
}




# predict.GAM_bamlss <- function(object, newdata, ...) {
#   x <- getS3method("predict", "bamlss")(object, newdata, ...)
#   return(x$mu)
# }




# Function from Stauffer
# getting the expectation of the model with transformation
powexp <- function(
    mean, sd, power = 1, left = 0, n = 5000, FUN = "mean"
) {

    if (!is.logical(left) && !missing(power)) {
        if (left < 0) {
            stop(
                sprintf(
                    paste0(
                        "Cannot easily apply power",
                        "parameter %f if left is lower than 0!"
                    ),
                    power
                )
            )
        }
    }
    applyfun <- function(x, resfun, p, left) {
        power <-  x[3L]
        x <- qnorm(p, x[1L], x[2L])
        i <- which(x > 0)
        x[i] <- x[i]^power
        if (!is.logical(left))
            x[which(x < left)] <- left
        if (is.character(resfun)) {
            if (resfun == "x") return(x)
        }
        do.call(resfun, list(x = x))
    }
    res <- data.frame(
        "mean" = mean, "sd" = sd, "power" = power, "exp" = rep(NA, length(mean))
    )
    p   <- seq(0, 1, length.out = (n + 2))[2:(n + 1)]
    apply(res, 1, applyfun, resfun = FUN, p = p, left = left)
}




# predict.GAM_bamlss_cnorm <- function(object, newdata, power=2, ...){
#   x <- getS3method("predict", "bamlss")(object, newdata, ...)
#   return(powexp(x$mu, exp(x$sigma), power = power))
# }