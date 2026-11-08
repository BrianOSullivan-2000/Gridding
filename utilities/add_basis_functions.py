
################################################################################

## Create basis functions at several resolutions as predictors for an ML      ##
## model. Generation uses a prespecified mesh (nodes) for spatial basis       ##
## functions, and time points for temporal basis functions.                   ##

## Can currently use Exponential, Matérn, or Wendland basis functions         ##
## Will also update a base formula with the basis functions                   ##

################################################################################

import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import gamma, kv


def add_basis_functions(
        coords, nodelist, supportlist = None,
        bfs_type = "Exp", parlist = (1, 50000, 0)
):
    
    # coords                -easting/northing coordinates of points
    # nodelist              -tuple of node coordinates for each mesh
    # supportlist           -tuple of support parameters
    # bfs_type              -type of basis function (Exp, Mat, Wendland)
    # parlist               -tuple of basis function parameters

    if supportlist is None:
        supportlist = [np.inf]
    basis_functions = np.zeros((coords.shape[0],
                                np.vstack(nodelist).shape[0]))
    
    for i in range(len(nodelist)):
        
        nodes = nodelist[i]
        sd = cdist(coords, nodes)
        bfs_pars = parlist[i]
        support = supportlist[i]

        if bfs_type == "Exp":
            sig2, phi, nugget = bfs_pars[0], bfs_pars[1], bfs_pars[2]
            res_basis_functions = (sig2 * np.exp(-sd / phi)) + nugget
        elif bfs_type == "Mat":
            sig2, phi, nugget, kappa = bfs_pars[0], bfs_pars[1], bfs_pars[2], bfs_pars[3]
            res_basis_functions = sig2 * ((2**(kappa - 1) * gamma(kappa))**-1) * ((sd / phi)**kappa) * kv(kappa, sd / phi) + nugget
        elif bfs_type == "Wendland":
            sig2, phi, nugget, power = bfs_pars[0], bfs_pars[1], bfs_pars[2], bfs_pars[3]
            res_basis_functions = np.zeros((sd.shape[0], sd.shape[1]))
            for p in range(power + 1):
                res_basis_functions += ((sd / phi)**p)
            res_basis_functions = sig2 * res_basis_functions + nugget
        elif bfs_type == "None":
            res_basis_functions = sd
            phi = 1

        res_basis_functions[(sd / phi) > support]  = 0

        start = sum(len(n) for n in nodelist[:i])
        end = sum(len(n) for n in nodelist[:i + 1])
        basis_functions[:, start:end] = res_basis_functions
        
    return basis_functions