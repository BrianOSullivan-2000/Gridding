
################################################################################

# Inverse Distance Weighting
# Simple interpolation method. Uses weighted average based on inverse power law

################################################################################

IDW <- function(y, coords, new_coords, idp = 2, nmax = 8){
  
  # y                       -observed values of response
  # coords                  -spatial coordinates of observations
  # new_coords              -spatial coordinates to interpolate on to
  #                         coords and new coords should have two columns
  #                         eg (east, north)
  
  # idp                     -inverse distance weighting power
  # nmax                    -maximum number of neighbours to use for prediction
  
  # Load in packages
  require(sp); require(gstat)
  
  coords <- as.matrix(coords); new_coords <- as.matrix(new_coords)
  
  # For gstat - the data needs to be stored in spatial dataframes
  data <- data.frame(y = y, east = coords[, 1], north = coords[, 2])
  coordinates(data) <- c("east", "north")
 
  newdata <- data.frame("east" = new_coords[, 1], "north" = new_coords[, 2])
  coordinates(newdata) <- c("east", "north")
  
  # Interpolate value at each point in newdata
  idw_df <- idw(y ~ 1, data, newdata, idp = idp, nmax = nmax,
                debug.level = 0)
  
  return(idw_df)
}

