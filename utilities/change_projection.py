import numpy as np
from pyproj import Transformer


def change_projection(coords, source_crs="EPSG:4326", target_crs="EPSG::29903"):

    transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)

    # Handle a single point
    if isinstance(coords, tuple):
        x, y = transformer.transform(coords[0], coords[1])
        return x, y

    # Handle an array of points
    coords = np.asarray(coords)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must have shape (N, 2)")

    x_new, y_new = transformer.transform(coords[:, 0], coords[:, 1])
    return np.column_stack((x_new, y_new))