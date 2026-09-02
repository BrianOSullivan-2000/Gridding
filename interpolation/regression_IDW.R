
################################################################################

# Inverse Distance Weighting on residuals from a regression method
# Simple interpolation method. Uses weighted average based on inverse power law.
# key step is you have to provide a regression formula

################################################################################

regression_IDW <- function(df, new_df, coords, new_coords,
                           regression_method,
                           idp = 2, nmax = 5, ...){
  
  # df                      -data frame of response variable and covariates
  # new_df                  -data frame of covariates for grid locations
  # coords                  -spatial coordinates of observations
  # new_coords              -spatial coordinates to interpolate on to
  #                         coords and new coords should have two columns
  #                         eg (east, north)
  
  # regression_method       -method to use for regression
  
  # idp                     -inverse distance weighting power
  # nmax                    -maximum number of neighbours to use for prediction
  
  # Load in packages
  require(sp); require(gstat)
  
  coords <- as.matrix(coords); new_coords <- as.matrix(new_coords)
  
  # Get the residuals
  model <- regression_method(df, coords = coords, ...)
  df$residuals <- model$residuals
  
  # For gstat - the data needs to be stored in spatial dataframes
  data <- data.frame(residuals = df$residuals, 
                     east = coords[, 1], north = coords[, 2])
  coordinates(data) <- c("east", "north")
  
  newdata <- data.frame("east" = new_coords[, 1], "north" = new_coords[, 2])
  coordinates(newdata) <- c("east", "north")

  # Interpolate value at each point in newdata
  idw_df <- idw(residuals ~ 1, data, newdata, idp = idp, nmax = nmax, debug.level = 0)
  new_df$residuals <- idw_df$var1.pred
  
  # Add trend back
  new_df$var1.pred <- predict(model, new_df) + idw_df$var1.pred
    
  return(new_df)
}

