
################################################################################

# Regression interpolation - just estimate gridded values from regression model

################################################################################

regression_grid <- function(df, new_df, regression_method, formula, ...){
  
  # df                      -data frame of response variable and covariates
  # new_df                  -data frame of covariates for grid locations
  
  # regression_method       -method to use for regression
  # formula                 -regression formula - will depend on method used
  
  # Fit the regression model
  model <- regression_method(df, formula = formula, ...)
  
  # Predict at grid points
  new_df$pred <- predict(model, new_df)
  
  return(new_df)
}


