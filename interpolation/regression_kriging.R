
################################################################################

# Ordinary Kriging
# Kriging interpolation - the response variable is considered to come from
# a Gaussian Process. Specify a covariance function (eg Matérn) and fit to 
# a variogram. Then interpolate using weighted averages.

################################################################################
  
regression_kriging <- function(df, new_df, coords, new_coords,
                               regression_method,
                               nmax = Inf, debug.level = 1, ...){
  
  # df                      -data frame of response variable and covariates
  # new_df                  -data frame of covariates for grid locations
  # coords                  -spatial coordinates of observations
  # new_coords              -spatial coordinates to interpolate on to
  #                         coords and new coords should have two columns
  #                         eg (east, north)
  
  # nmax                    -maximum number of neighbours to use for prediction
  # debug.level             -set to -1 to see progress of kriging
  #                         set to 0 for no printed output
  
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
  
  # Fit spatial variogram to values
  fit_vgm <- spatial_variogram(data, ...)
  
  # Interpolate values with ordinary kriging
  krig_df <- krige(residuals ~ 1, data, newdata, fit_vgm, nmax = nmax, 
                   debug.level = debug.level)
  new_df$residuals <- krig_df$var1.pred
  
  # Add trend back
  new_df$var1.pred <- predict(model, new_df) + krig_df$var1.pred
  
  return(new_df)
}




## Make and fit a spatial variogram

spatial_variogram <- function(df, flex_vgm = F, flex_fit = T,
                              cutoff = NA, width = NA,
                              vgm_model = "Mat",
                              psill = NA, nugget = NA, range = NA, 
                              kappa = 1, fit.method = 7, plot_vgm = F, ...){
  
  # df               -data
  # flex_vgm         -automate making empirical variogram (Not recommended)
  # flex_fit         -automate initial values for theoretical variogram
  #                  (Not recommended)
  
  # cutoff           -cutoff of empirical variogram
  # width            -width of bins for empirical variogram
  
  # vgm_model        -variogram model
  # psill, nugget,
  # range            -initial partial sill, nugget and range
  # kappa            -kappa value (shape parameter for Matern model)
  # fit.method       -fitting method for vgm (see fit.variogram() documentation)
  # plot_vgm         -option to plot the variogram
  
  # Load in packages
  require(sp); require(gstat)
  
  # Error message if variogram params not provided and flex arguments are False
  if ((any(is.na(c(cutoff, width))) & flex_vgm == F) |
      (any(is.na(c(psill, nugget, range))) & flex_fit == F)){
    
    # Find parameters that are missing
    params <- c("cutoff", "width",
                "psill", "nugget", "range")[is.na(c(cutoff, width,
                                                    psill, nugget, range))]
    params <- paste(params, collapse = ", ")
    
    warning <- paste0("No values assigned to variogram parameters. 
    Either assign parameters manually (recommended) or use flex_vgm and
    flex_fit arguments. Parameters missing are ", params, ".")
    
    stop(warning)
  }
  
  # Flexible way to get cutoff and width of variogram
  if (flex_vgm == T){
    cutoff <- sqrt((max(df$east) - min(df$east))^2 + 
                     (max(df$north) - min(df$north))^2) / 2
    width = cutoff / 15
  }
  
  # Empirical variogram
  vgm <- variogram(residuals ~ 1, data = df, cutoff = cutoff, width = width)
  
  # Flexible way to get initial variogram parameters 
  if (flex_fit == T){
    nugget <- mean(vgm$gamma[1:3])
    psill <- mean(vgm$gamma) - nugget
    range <- vgm$dist[which.min(abs(vgm$gamma - (psill + nugget)))]
  }
  
  # Set up theoretical variogram - then fit to empirical variogram
  fit_vgm <- vgm(model = vgm_model, psill = psill, nugget = nugget,
                 range = range, kappa = kappa)
  fit_vgm <- fit.variogram(vgm, fit_vgm, fit.method = fit.method)
  
  if(plot_vgm){
    print(plot(vgm, fit_vgm))
  }
  
  return(fit_vgm)
}
