
################################################################################

# Spatio-temporal Inverse Distance Weighting
# Just like normal inverse distance weighting, except "temporal distance"
# is also considered. Interpolation done on residuals after regression

################################################################################

spatio_temporal_regression_IDW <- function(df, new_df, coords, new_coords,
                                           regression_method, 
                                           C = 500000, idp = 2, nmax = 8, ...){
  
  # df                      -data frame of response variable and covariates
  # new_df                  -data frame of covariates for grid locations
  # coords                  -coordinates of observations (east, north, time)
  # new_coords              -coordinates to interpolate on to (east, north, time)
  #                         coords and new coords should have two columns
  #                         eg (east, north)
  
  # regression_method       -method to use for regression
  # C                       -anisotropy, relates spatial and temporal distance
  #                         if default, just do spatial IDW
  # idp                     -inverse distance weighting power
  # nmax                    -maximum number of neighbours to use for prediction
  
  # Load in packages
  require(sp); require(gstat)
  
  coords <- as.matrix(coords); new_coords <- as.matrix(new_coords)
  
  # Get the residuals
  model <- regression_method(df, coords = coords, ...)
  df$residuals <- model$residuals
  
  # For gstat - the data needs to be stored in spatial dataframes
  data <- data.frame(residuals = df$residuals, east = coords[, 1],
                     north = coords[, 2], t = coords[, 3])
  newdata <- data.frame("east" = new_coords[, 1], 
                        "north" = new_coords[, 2], "t" = new_coords[, 3])
  
  data$t <- data$t * C; newdata$t <- newdata$t * C
  coordinates(data) <- c("east", "north", "t")
  coordinates(newdata) <- c("east", "north", "t")
  
  # Interpolate value at each point in newdata
  idw_df <- idw(residuals ~ 1, data, newdata, idp = idp, nmax = nmax)
  new_df$residuals <- idw_df$var1.pred
  
  # Add trend back
  new_df$var1.pred <- predict(model, new_df) + idw_df$var1.pred  

  return(new_df)
}
