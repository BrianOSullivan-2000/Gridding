
################################################################################

# fdaPDE - function data analysis method

# Developed by Laura Sangalli - the best reference for this is from 2021
# https://doi.org/10.1111/insr.12444

# Models data as a piecewise sum of linear basis functions over a mesh
# The power of fdaPDE is using a PDE that regularises the fit
# However, this function is the simplest implementation of it without any PDE

# You can install this on their GitHub: https://github.com/fdaPDE
# But BE WARNED
# Installing this thing was a nightmare, the developers are constantly updating
# their core repository but they haven't updated their C++ or R
# repositories to comply with it. (at least as of September 4th 2025)
# The normal install (ie devtools) doesn't work for this reason

# There is an older version of the package fdaPDE on CRAN too, I have little
# experience with this but it might be worth trying as an alternative

# To save yourself some time (it took me a few days to figure this out)
# just contact me at Brian.OSullivan@met.ie - I can help out

################################################################################

fdaPDE_noPDE <- function(df, new_df, coords, new_coords,
                         outline, mesh, formula,
                         old_crs = 29903, new_crs = 4326, bs_type = "P1",
                         a = 0.02, q = 30, lambda = NA,
                         lambda_grid = 10^seq(from = -8, to = -4, by = 0.5)) {

    # df                      -data frame of response variable and covariates
    # new_df                  -data frame of covariates for grid locations
    # coords                  -spatial coordinates of observations
    # new_coords              -spatial coordinates to interpolate on to
    #                         coords and new coords should have two columns
    #                         eg (east, north)

    # outline                 -outside boundary for the mesh
    # mesh                    -mesh object with vertices for basis functions
    # bs_type                 -type of basis function to use (polynomial 1 or 2)
    # a                       -maximum area allowed for triangles in mesh
    # q                       -minimum angle for triangles allowed in mesh

    # lambda                  -smoothness parameter - if NA will use gcv
    # lambda_grid             -candidate lambda values to try in gcv

    # Load in packages
    require(fdaPDE2)
    require(sf)
    require(RTriangle)
    require(mapview)
    require(dplyr)

    ## Before doing anything - need to convert easting/northing to lon/lat
    coords <- st_as_sf(coords, coords = names(coords), crs = old_crs)
    coords <- st_transform(coords, new_crs)
    df[c("lon", "lat")] <- st_coordinates(coords)

    new_coords <- st_as_sf(
        new_coords,
        coords = names(new_coords),
        crs = old_crs
    )
    new_coords <- st_transform(new_coords, new_crs)
    new_df[c("lon", "lat")] <- st_coordinates(new_coords)

    # If no mesh provided you need to make one
    if (missing(mesh)) {

        # Load in or create a boundary
        if (missing(outline)) {
            # If no outline given, just make a bounding box around the data
            outline <- st_bbox(coords)
            pad_dist <-
                0.01 * max(
                    (outline["xmax"] - outline["xmin"]),
                    (outline["ymax"] - outline["ymin"])
                )
            outline <- st_buffer(st_as_sfc(outline), dist = pad_dist)
        } else if (simplify_outline) {
            outline <- concaveman(outline, concavity = 10)
            pad_dist <- 0.01 * max(
                (outline["xmax"] - outline["xmin"]),
                (outline["ymax"] - outline["ymin"])
            )
            outline <- st_buffer(outline, dist = pad_dist)
            outline <- st_simplify(outline, dTolerance = pad_dist)
        } else if (nrow(st_coordinates(outline)) > 60) {
            warning(
                paste0(
                    "The boundary provided is quite complex (> 60 vertices). ",
                    "Consider using concaveman() and st_simplify() ",
                    "to reduce its complexity."
                )
            )
        }

        if (q < 30) {
            warning(
                paste0(
                    "Acute angles in the mesh may lead to ",
                    "numerical instability. Consider q >= 30"
                )
            )
        }

        # Get coordinates of outline vertices
        outline_vertices <- st_cast(x = outline, "POINT", crs = 4326)
        outline_vertices <- st_coordinates(x = outline_vertices)
        outline_vertices <- data.frame(
            lon = outline_vertices[, 1],
            lat = outline_vertices[, 2]
        )
        outline_vertices <- outline_vertices[-nrow(outline_vertices), ]

        # Define segments
        outline_segments <- cbind(
            seq_len(nrwo(outline_vertices)), c(2:nrow(outline_vertices), 1)
        )

        # Make the mesh
        outline_pslg <- pslg(P = outline_vertices, S = outline_segments)
        mesh <- triangulate(p = outline_pslg, a = a, q = q)
        mesh$H <- matrix(data = numeric(0), ncol = 2)
        mesh <- triangulation(
            nodes = mesh$P,
            cells = mesh$T,
            boundary = mesh$PB
        )
    }

    # Make a geoframe and add data as a layer
    gf <- geoframe(domain = mesh)
    df_vars <- c(all.vars(formula), "lon", "lat")
    df <- as.data.frame(df |> dplyr::select(df_vars))

    gf$insert(
        layer = "response", type = "point",
        geo = c("lon", "lat"), data = df
    )

    ## Set up basis functions - P1 or P2
    basis_functions <- fe_function(domain = mesh, type = bs_type)

    ## You have to add basis_functions to the global environment due to how
    ## the package is coded up. If you don't like this, fair, just contact me
    ## and we can see if there is an alternative way to fix it
    assign("basis_functions", basis_functions, envir = globalenv())

    ## Add basis functions to formula
    fda_formula <- update(formula, . ~ . + basis_functions)

    ## Set up spatial regression model note the formula needs to be updated to
    model <- sr(formula = fda_formula, data = gf)

    ## Fit model (either fit fixed lambda or gcv grid search)
    if (is.na(lambda)) {
        model_fit <-
            model$fit(
                calibrator = gcv(optimizer = grid_search(grid = lambda_grid))
            )
        print(model_fit$optimum)
    } else {
        model_fit <- model$fit(lambda = lambda)
    }

    ## Evaluate model at new locations
    new_df["pred"] <- basis_functions$eval(new_df[c("lon", "lat")])

    ## Add covariate effect if included
    if (length(all.vars(formula)) > 1) {
        new_df["pred"] <-
            new_df["pred"] +
            (as.matrix(
                new_df[attr(terms(formula), "term.labels")]
            ) %*% model$beta)
    }

    new_df
}
