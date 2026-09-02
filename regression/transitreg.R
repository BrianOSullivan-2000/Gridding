
################################################################################

# Transitional Regression Model (transitreg)
# Developed by Reto Stauffer (2025). By binning values, one can fit
# continuous data using a transitional regression model 

################################################################################

TransitReg <- function(df, formula, breaks,
                       plot = F, readout = F, ...){
  
  # df                      -data frame of response variable and covariates
  # formula                 -regression formula
  # breaks                  -define breaks for binning continuous data
  
  # plot                    -plot the output of the regression model
  # readout                 -print a summary of the model
  
  # transitreg package is available on GitHub
  require("transitreg")
  
  model <- transitreg(formula = formula, data = df, breaks = breaks, ...)
  class(model) <- c("TransitReg", class(model)) 
  
  # Plot and summary
  if (plot == T){
    plot(model)
  }
  if (readout == T){
    print(summary(model))
  }
  
  return(model)
}




predict.TransitReg <- function(object, newdata, ...){
  x <- getS3method("predict", "transitreg")(object, newdata, 
                                            type = "quantile", p = 0.5, ...)
  return(x)
}