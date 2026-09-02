
################################################################################

# Multiple linear regression function

################################################################################

MLR <- function(df, formula, step_method = "forward",
                plot = F, readout = F, ...){
  
  # df                      -data frame of response variable and covariates
  # formula                 -regression formula
  # step_method             -method to carry out stepwise regression
  #                         (forward, backward, both, none)
  #                         if none, no step-wise regression will be calculated
  
  # plot                    -plot the output of the regression model
  # readout                 -print a summary of the model
  
  # MASS package required for step-wise
  require(MASS)
  
  # Standard Regression  
  model <- lm(formula = formula, data = df)
  
  if(tolower(step_method) != "none"){
    
    # Stepwise Regression
    model <- stepAIC(model, method = step_method, trace = F)
  }
  
  # Plot and Summary
  if(plot == T){plot(model)}
  if(readout == T){print(summary(model))}
  
  return(model)
}