################################################################################

## Add basis functions at several resolutions as predictors to a dataframe    ##
## using a prespecified mesh, mostly for seRF and DeepKriging                 ##

## Can currently use Exponential, Matérn, or Wendland basis functions         ##
## Will also update a base formula with the basis functions                   ##

################################################################################

add_basis_functions <- function(
        coords, nodelist, supportlist = list(Inf),
        bfs_type = "Exp", parlist = list(c(1, 20, 0)),
        base_formula = "y ~ 1",
        center_params = attributes(scale(coords))$`scaled:center`,
        scale_params = attributes(scale(coords))$`scaled:scale`) {

    # coords                -easting/northing coordinates of points
    # nodelist              -list of nodes at each resolution
    # supportlist           -list of support parameters

    # bfs_type              -type of basis function (Exp, Mat, Wendland)
    # parlist               -list of basis function parameters
    # base_formula          -regression formula before adding basis functions
    # center_params         -center parameters for scaling coordinates
    # scale_params          -scale parameters for scaling coordinates

    require(sp)

    basis_functions <- matrix(0,
        nrow = nrow(coords),
        ncol = sum(unlist(lapply(nodelist, nrow)))
    )

    ## Loop through each resolution
    for (i in seq_along(nodelist)) {

        nodes <- nodelist[[i]]
        coords_normalized <- coords

        ## Normalize the nodes and coordinates
        nodes[, 1] <- (nodes[, 1] - center_params[1]) / scale_params[1]
        nodes[, 2] <- (nodes[, 2] - center_params[2]) / scale_params[2]
        coords_normalized[, 1] <-
            (coords[, 1] - center_params[1]) / scale_params[1]
        coords_normalized[, 2] <-
            (coords[, 2] - center_params[2]) / scale_params[2]

        ## Distances, parameters, and support
        sd <- spDists(coords_normalized, nodes)
        bfs_pars <- parlist[[i]]
        support <- supportlist[[i]]

        ## Depending on bfs_type, make the basis functions
        if (bfs_type == "Exp") {

            sig2 <- bfs_pars[1]
            phi <- bfs_pars[2]
            nugget <- bfs_pars[3]
            res_basis_functions <- (sig2 * exp(-sd / phi)) + nugget

        } else if (bfs_type == "Mat") {

            require(geoR)
            sig2 <- bfs_pars[1]
            phi <- bfs_pars[2]
            nugget <- bfs_pars[3]
            kappa <- bfs_pars[4]
            res_basis_functions <-
                sig2 *
                ((2^(kappa - 1) * gamma(kappa))^(-1)) *
                ((sd / phi)^kappa) *
                besselK(x = sd / phi, nu = kappa) + nugget

        } else if (bfs_type == "Wendland") {

            sig2 <- bfs_pars[1]
            phi <- bfs_pars[2]
            nugget <- bfs_pars[3]
            power <- bfs_pars[4]
            res_basis_functions <- matrix(0, nrow = nrow(sd), ncol = ncol(sd))
            for (p in 0:power) {
                res_basis_functions <- res_basis_functions + ((sd / phi)^p)
            }
            res_basis_functions <- sig2 * res_basis_functions + nugget

        } else if (bfs_type == "None") {
            res_basis_functions <- sd
            phi <- 1
        }

        ## Remove any values that are beyond the support
        res_basis_functions[(sd / phi) > support] <- 0

        ## Update basis functions with ones from current resolution
        basis_functions[, (1 + sum(unlist(lapply(nodelist[0:(i - 1)], nrow)))):
                            sum(unlist(lapply(nodelist[0:i], nrow)))] <-
            res_basis_functions
    }

    ## Add column names and add basis functions to formula
    colnames(basis_functions) <- paste0("B", seq_along(basis_functions))
    f <- as.formula(paste(base_formula, "+",
                          paste(colnames(basis_functions), collapse = "+")))

    list("basis_functions" = basis_functions, "f" = f)
}
