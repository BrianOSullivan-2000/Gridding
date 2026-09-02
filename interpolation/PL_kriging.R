################################################################################

# Kriging with a model fit by Penalised Likelihood
# The model is still just normal ordinary/universal kriging
# However, the model is estimated using Penalised Likelihood. We're adding
# an Elastic-net penalty (via glmnet) but others could be considered in future

# The model is fit in a two-step process. First the covariance structure is fit,
# then the regression parameters (beta) are fit with penalized regression. We 
# iterate between these two steps, then do kriging with the final model

# In addition - I've included the option to constrain regression parameters

################################################################################

PL_kriging <- function(df, new_df, coords, new_coords,
                       
                       init_pars, lower, upper,
                       formula = y ~ 1, cov_function = "Mat", maxit = 1000,
                       
                       lambda = 0.1, alpha = 0.5, 
                       constrained, beta_lower, beta_upper,
                       
                       max_loops = 3, perc_tol = 1e-8,
                       
                       useRcpp = F){
  
  # df                      -data frame of response variable and covariates
  # new_df                  -data frame of covariates for grid locations
  # coords                  -coordinates of observations (east, north)
  # new_coords              -coordinates to interpolate on to (east, north)
  #                         coords and new coords should have two columns
  #                         eg (east, north)
  
  # init_pars               -initial covariance params
  # lower/upper             -lower and upper limits of covariance params
  # formula                 -regression formula (linear regression)
  # cov_function            -type of covariance function to use, can use
  #                         Exponential("Exp"), Matérn("Mat"), Spherical("Sph")
  # maxit                   -max number of iterations for MLE
  
  # lambda                  -scaling parameter for penalty (how strong is it)
  # alpha                   -describes how much of the regression is lasso 
  #                         (alpha = 1) vs. ridge (alpha = 0)
  # constrained             -logic vector, label constrained regression params
  # beta_lower/beta_upper   -set limits for constrained params
  
  # max_loops               -the max number of iterations between fitting 
  #                         the covariance and fitting regression params
  # perc_tol                -percentage tolerance between params for convergence
  
  # useRcpp                 -if true, will use faster functions coded in C++
  #                         (Just for likelihood + kriging, still uses glmnet)
  
  # Load in necessary packages
  require(sp); require(dplyr); require(optimParallel)
  
  ## Setup for optimized code, set up cluster for parallelization, 
  ## and source all the Rcpp compiled C++ code on that cluster
  if(useRcpp){solver <- optimParallel}
  else{solver <- optim}
  solver <- optim
  
  
  ## Create spatial objects and get distances between points
  sdf <- df; new_sdf <- new_df
  sdf$east <- coords[, 1]; sdf$north <- coords[, 2]
  new_sdf$east <- new_coords[, 1]; new_sdf$north <- new_coords[, 2]
  coordinates(sdf) <- c("east", "north")
  coordinates(new_sdf) <- c("east", "north")
  ds <- spDists(sdf)
  cross_ds <- spDists(sdf, new_sdf)
  
  
  ## Set up parameters for covariance function
  if(cov_function %in% c("Exp", "Gau", "Sph")){
    if(missing(init_pars)){
      pars = c(1.5, 0.5, 50000)
      lower = c(0.001, 0, 1000)
      upper = c(Inf, Inf, Inf) 
    }
    else if (length(init_pars) != 3){stop("Wrong number of parameters (should be three)")}
  }
  if(cov_function == "Mat"){
    if(missing(init_pars)){
      pars = c(1.5, 0.5, 50000, 1.5)
      lower = c(0.001, 0, 1000, 0.3)
      upper = c(Inf, Inf, Inf, 2.5)
    }
    else if (length(init_pars) != 4){stop("Wrong number of parameters (should be four)")}
  }
  if(!(cov_function %in% c("Exp", "Gau", "Sph", "Mat"))){
    stop("Model is not available, choose either Exp, Gau, Sph, or Mat")
  }
  if(!missing(init_pars)){pars <- init_pars}
  
  
  ## First pass - fit the covariance parameters using normal MLE
  ## Response variable and covariates
  y <- model.frame(formula, df)[[1]]; X <- model.matrix(formula, df)
  new_X <- model.matrix(reformulate(attr(terms(formula), 
                                         "term.labels")), new_df)

  ## Optimize likelihood numerically to get covariance parameters
  soln <- solver(par = pars, fn = spatial_likelihood, ds = ds,
                 y = y, X = X, cov_function = cov_function,
                 lower = lower, upper = upper, 
                 control = list(maxit = maxit, parscale = pmax(pars, 1e-4), factr = 1e8),
                 method = "L-BFGS-B", useRcpp = useRcpp)
  pars <- soln$par
  
  ## Initial setup for iterative loop
  perc_diff <- rep(Inf, length(pars))
  current_loop <- 0
  
  ## Iteration ends either at max loops or at convergence
  while((current_loop < max_loops) & (max(perc_diff) > perc_tol)){
    
    ############################################################################
    ## First part of loop - regression                                        ##
    ############################################################################
    
    ## Get current covariance matrix
    if(useRcpp){
      ## For Rcpp covariance functions, need to manipulate the arguments a little
      fn <- get(paste0("cov_", cov_function, "_cpp"))
      names(pars) <- names(formals(func))[2:(1 + length(pars))]
      V <- do.call(fn, c(list(ds = ds), as.list(pars), list(diagonal = T)))
    }
    else{
      fn <- get(paste0("cov_", cov_function))
      V <- fn(pars, ds)
    }
    
    ## Scale y and X
    R <- chol(V)
    weighted_y <- backsolve(R, y); weighted_X <- backsolve(R, X)
    
    ## Fit beta parameters using Elastic-Net regularisation
    model <- glmnet(weighted_X[, -1], weighted_y, standardize = F,
                    lambda = lambda, alpha = alpha)
    beta <- c(as.numeric(model$a0), as.numeric(model$beta))
          
    ## Adjusting according to constraints
    if(!missing(constrained)){
      constrained_idx <- which(constrained)
      
      ## Update betas that are outside limits
      previous_beta <- beta
      beta[constrained_idx] <- 
        pmin(pmax(beta[constrained_idx], beta_lower), beta_upper)
      
      ## Drop columns in X corresponding to fixed beta values
      changed_idx <- which(beta != previous_beta)
      if(length(changed_idx) > 0){
        
        need_beta_update <- T
        reduced_X <- weighted_X[, -changed_idx, drop = FALSE]
        
        ## Redo the elastic-net with new X, update the unconstrained beta params
        model <- glmnet(reduced_X[, -1], weighted_y, standardize = F,
                        lambda = lambda, alpha = alpha)
        beta[setdiff(seq_along(beta), changed_idx)] <- model$beta
      }
      
      while(need_beta_update){
        
        ## Check all betas again
        previous_beta <- beta
        beta[constrained_idx] <- 
          pmin(pmax(beta[constrained_idx], beta_lower), beta_upper)
        
        if(sum(beta != previous_beta) > 0){
          
          changed_idx <- c(changed_idx, which(beta != previous_beta))
          reduced_X <- weighted_X[, -changed_idx, drop = FALSE]
          
          ## Redo the elastic-net with new X, update the unconstrained beta params
          model <- glmnet(reduced_X[, -1], weighted_y, standardize = F,
                          lambda = lambda, alpha = alpha)
          beta[setdiff(seq_along(beta), changed_idx)] <- model$beta
        }
        else{
          need_beta_update <- F
        }
      }
    }
    
    ############################################################################
    ## Second part of loop - fit covariance (Likelihood)                      ##
    ############################################################################
    soln <- solver(par = pars, fn = spatial_likelihood, ds = ds,
                   y = y, X = X, cov_function = cov_function,
                   lower = lower, upper = upper, 
                   control = list(maxit = maxit, parscale = pmax(pars, 1e-4), factr = 1e8),
                   method = "L-BFGS-B", beta = beta, useRcpp = useRcpp)
    previous_pars <- pars
    pars <- soln$par
    
    perc_diff <- abs((pars - previous_pars)/pmax(previous_pars, 1e-12))
    current_loop <- current_loop + 1
  }

  
  ## Final Step - interpolate using kriging
  if(useRcpp){
    if(is.null(beta)){ beta <- numeric(0)}
    new_df$pred <- predict_MVN_cpp(pars = pars, ds = ds, cross_ds = cross_ds,
                                   y = y, X = X, new_X = new_X,
                                   cov_function = cov_function, beta = beta)
  }
  else{
    # Get covariance function, inverse, and determinant
    fn <- get(paste0("cov_", cov_function))
    V <- fn(pars, ds)
    cross_V <- fn(pars, cross_ds)
    Vinv <- solve(V)
    
    new_df$pred <- as.numeric(new_X %*% beta + t(cross_V) %*%
                                (Vinv %*% (y - X %*% beta)))
  }
  
  print(cbind(terms, beta))
  return(new_df)
}




## For the record, this function is identical to spatial_likelihood in MLE_kriging.R
spatial_likelihood <- function(pars, ds, y, X, cov_function, 
                               beta = NULL, useRcpp = F){
  
  if(useRcpp){
    if(is.null(beta)){ beta <- numeric(0)}
    negloglik <- spatial_likelihood_cpp(pars = pars, ds = ds, y = y, X = X,
                                        cov_function = cov_function, beta = beta)
  }
  
  else{
    # Get covariance function, inverse, and determinant
    fn <- get(paste0("cov_", cov_function))
    V <- fn(pars, ds)
    
    Vchol <- chol(V); logdetV <- 2 * sum(log(diag(Vchol)))
    
    # Get regression parameters
    if(missing(beta)){
      Vinvy <- backsolve(Vchol, forwardsolve(t(Vchol), y))
      VinvX <- backsolve(Vchol, forwardsolve(t(Vchol), X))
      
      betahat <- solve(t(X) %*% VinvX) %*% (t(X) %*% Vinvy)
    }
    else{betahat <- beta}
    
    res <- y - (X %*% betahat)
    Vinvres <- backsolve(Vchol, forwardsolve(t(Vchol), res))
    
    # Calculate negative log-likelihood
    negloglik <- 0.5 * (logdetV + (t(res) %*% Vinvres))
  }
  return(negloglik)
}




## The four different covariance functions available
## Again these are also in MLE_kriging.R
cov_Exp <- function(pars, ds){
  sig2 <- pars[1]; nugget <- pars[2]; phi <- pars[3]
  V <- sig2*exp(-ds/phi)
  V[ds == 0] <- V[ds == 0] + (sig2 + nugget)
  return(V)
}
cov_Sph <- function(pars, ds) {
  sig2   <- pars[1]; nugget <- pars[2]; phi <- pars[3]
  h <- ds/phi; idx <- h <= 1; V <- ds*0
  
  V[idx] <- sig2 * (1 - 1.5*h[idx] + 0.5*h[idx]^3)
  V[ds == 0] <- V[ds == 0] + (sig2 + nugget)
  return(V)
}
cov_Gau <- function(pars, ds){
  sig2 <- pars[1]; nugget <- pars[2]; phi <- pars[3]
  V <- sig2*exp(-(ds/phi)^2)
  V[ds == 0] <- V[ds == 0] + (sig2 + nugget)
  return(V)
}
cov_Mat <- function(pars, ds){
  sig2 <- pars[1]; nugget <- pars[2]; phi <- pars[3]; kappa = pars[4]
  h <- ds/phi; idx <- h > 0; V <- ds*0
  
  V[idx] <- sig2 * ((2^(kappa-1)*gamma(kappa))^(-1)) *
    ((h[idx])^kappa) *
    besselK(x = h[idx], nu = kappa)
  V[ds == 0] <- V[ds == 0] + (sig2 + nugget)
  return(V)
}
