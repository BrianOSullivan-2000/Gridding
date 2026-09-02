
################################################################################

# Spatio-temporal Inverse Distance Weighting
# Just like normal inverse distance weighting, except distance 
# over time is also considered

################################################################################

spatio_temporal_IDW <- function(y, coords, new_coords, 
                                C = 500000, idp = 2, nmax = 8){
  
  # y                       -observed values of response
  # coords                  -coordinates of observations (east, north, time)
  # new_coords              -coordinates to interpolate on to (east, north, time)
  #                         coords and new coords should have two columns
  #                         eg (east, north)
  
  # C                       -anisotropy, relates spatial and temporal distance
  #                         if default, just do spatial IDW
  # idp                     -inverse distance weighting power
  # nmax                    -maximum number of neighbours to use for prediction
  
  # Load in packages
  require(sp); require(gstat)
  
  coords <- as.matrix(coords); new_coords <- as.matrix(new_coords)
  
  # For gstat - the data needs to be stored in spatial dataframes
  data <- data.frame(y = y, east = coords[, 1],
                     north = coords[, 2], t = coords[, 3])
  newdata <- data.frame("east" = new_coords[, 1], 
                        "north" = new_coords[, 2], "t" = new_coords[, 3])
  
  data$t <- data$t * C; newdata$t <- newdata$t * C
  

  coordinates(data) <- c("east", "north", "t")
  coordinates(newdata) <- c("east", "north", "t")

  
  # Interpolate value at each point in newdata
  idw_df <- idw(y ~ 1, data, newdata, idp = idp, nmax = nmax)
  
  return(idw_df)
}
