################################################################################

# Random Forest Spatial Prediction Framework (RFsp)

# Standard ML algorithm using Random Forest. Spatial effects are embedded by
# including distances to all observations as covariates. Developed by Hengl.

# http://dx.doi.org/10.7717/peerj.5518

# This is like a primitive to the spatial embedded Random Forests (seRF)
# that I have developed

################################################################################

## Need the random forest and regression_grid functions
.file <- sys.frame(1)$ofile
.dir <- dirname(normalizePath(.file, mustWork = TRUE))
source(file.path(.dir, "..", "regression", "random_forest.R"))
source(file.path(.dir, "regression_grid.R"))
rm(.file, .dir)

RFsp <- function(
    df, new_df,
    coords, new_coords,
    base_formula = "y ~ 1",
    sd = NULL, new_sd = NULL,
    ...
) {
    # df                      -data frame of response variable and covariates
    # new_df                  -data frame of covariates for grid locations

    # coords                  -spatial coordinates of observations
    # new_coords              -spatial coordinates to interpolate on to
    #                         coords and new coords should have two columns
    #                         eg (east, north)

    # base_formula            -regression formula of all other covariates
    #                         (e.g. elevation)

    # sd, new_sd              -distance matrices, used to determine spatial
    #                         predictors. Normally these will be calculated,
    #                         but can be added in advance to save compute

    require(sp)
    require(dplyr)

    ## Estimates are more stable if scale to km
    coords <- coords |>
        mutate(east = east / 1000, north = north / 1000) |>
        as.matrix()
    new_coords <- new_coords |>
        mutate(east = east / 1000, north = north / 1000) |>
        as.matrix()

    ## Calculate distance matrices (for spatial predictors)
    if (is.null(sd) || is.null(new_sd)) {
        sd <- spDists(coords)
        new_sd <- spDists(new_coords, coords)
        colnames(sd) <- paste0("D", seq_len(ncol(sd)))
        colnames(new_sd) <- paste0("D", seq_len(ncol(new_sd)))
    }

    ## Add distances to data
    df <- cbind(df, sd)
    new_df <- cbind(new_df, new_sd)

    ## Update regression formula
    f <- as.formula(paste(
        base_formula, "+", paste(colnames(sd), collapse = "+")
    ))

    ## Return predicted values
    regression_grid(
        df,
        new_df,
        regression_method = Random_Forest,
        formula = f,
        ...
    )
}