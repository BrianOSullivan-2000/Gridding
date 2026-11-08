import matplotlib.pyplot as plt
import numpy as np
from meshpy.triangle import MeshInfo, build
from shapely import concave_hull
from shapely.geometry import MultiPoint


def make_mesh(island_outline, a, q, buffer=5000, simplify=10000, plot=False):
    
    ## Make boundary a simpler polygon
    points = MultiPoint(list(zip(island_outline['x'], island_outline['y'])))
    boundary = concave_hull(points, ratio=0.3)
    boundary = boundary.buffer(buffer)
    boundary = boundary.simplify(simplify)

    # Extract boundary coordinates and scale them
    coords = np.array(list(boundary.exterior.coords)[:-1])
    mean = coords.mean(axis=0)
    scale = coords.std(axis=0)
    scaled_coords = (coords - mean) / scale

    # Nodes and Edges (info for mesh)
    points = [(x, y) for x, y in scaled_coords]
    facets = [
        (i, (i + 1) % len(points))
        for i in range(len(points))
    ]
    mesh_info = MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(facets)

    # Generate mesh
    mesh = build(
        mesh_info,
        max_volume=a,
        min_angle=q
    )
    mesh_points = np.array(mesh.points)
    mesh_points = np.array(mesh.points) * scale + mean

    if plot:
        plt.figure()

        for element in mesh.elements:
            tri = mesh_points[element]
            plt.plot(
                [tri[0, 0], tri[1, 0], tri[2, 0], tri[0, 0]],
                [tri[0, 1], tri[1, 1], tri[2, 1], tri[0, 1]]
            )

        plt.scatter(island_outline["x"], island_outline["y"], c="grey", s=1)
        plt.axis("equal")
        plt.show()

    return(mesh_points)