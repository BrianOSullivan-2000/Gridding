
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator


def base_grid_plot(grid, stations=None,
                   response="y", bbox=None, 

                   background=None, background_stations=None,
                   
                   title="", title_size=10,
                   legend_title="", legend_title_size=10,

                   color_palette=None, bias = 1,
                   breaks=None, nbreaks=10, ymin=None, ymax=None,

                   station_size = 50,

                   fig_size=(12, 12), save_file=None, ax=None,

                   include_cbar=True, include_axis=True):

    ## Make the raster

    if color_palette is None:
        color_palette = ["#ffffd9", "#1d91c0", "#081d58"]
    
    east_coords = np.sort(grid["east"].unique())
    north_coords = np.sort(grid["north"].unique())

    raster = (
        grid.pivot(index="north", columns="east", values=response)
        .reindex(index=north_coords, columns=east_coords)
        .values
    )

    ## Origin starts at top left, so we have to flip the raster
    raster = np.flipud(raster)
    extent = (
        east_coords.min(),
        east_coords.max(),
        north_coords.min(),
        north_coords.max()
    )

    ## Set up color scheme

    if breaks is None:
        if ymin is None or ymax is None:
            ymin = grid[response].min(skipna=True)
            ymax = grid[response].max(skipna=True)

        locator = MaxNLocator(nbins=nbreaks)
        breaks = locator.tick_values(ymin, ymax)

    positions = np.linspace(0, 1, len(color_palette))**bias

    custom_cmap = mcolors.LinearSegmentedColormap.from_list(
        "custom_palette",
        list(zip(positions, color_palette))
    )

    color_palette = [
        mcolors.to_hex(custom_cmap(x))
        for x in np.linspace(0, 1, len(breaks) - 1)
    ]

    cmap = mcolors.ListedColormap(color_palette)
    norm = mcolors.BoundaryNorm(
        breaks,
        ncolors=len(color_palette),
        clip=True
    )

    new_plot = False
    if ax is None:
        _fig, ax = plt.subplots(figsize=fig_size)
        new_plot = True

    ## Can add a background
    if background is not None:

        bg_east = np.sort(background["east"].unique())
        bg_north = np.sort(background["north"].unique())

        bg_extent = (
            bg_east.min(),
            bg_east.max(),
            bg_north.min(),
            bg_north.max()
        )

        bg_raster = (
            background.assign(_bg=1)
            .pivot(index="north", columns="east", values="_bg")
            .reindex(index=bg_north, columns=bg_east)
            .values
        )

        bg_raster = np.flipud(bg_raster)


        ax.imshow(
            bg_raster,
            extent=bg_extent,
            cmap=mcolors.ListedColormap(["lightgrey"]),
            zorder=0
        )

    if background_stations is not None:
        ax.scatter(
            background_stations["east"],
            background_stations["north"],
            s=50,
            c="lightgrey",
            edgecolors="black",
            linewidth=0.5,
            alpha=0.5,
            zorder=1
        )


    img = ax.imshow(
        raster,
        extent=extent,
        cmap=cmap,
        norm=norm
    )

    if stations is not None:
        ax.scatter(
            stations["east"],
            stations["north"],
            c=stations[response],
            cmap=cmap,
            norm=norm,
            s=station_size,
            edgecolors="black",
            linewidth=0.5,
            zorder=10
        )

    if include_cbar:
        cbar = plt.colorbar(img, ax=ax, ticks=breaks)
        cbar.set_label(legend_title, fontsize=legend_title_size)

    if bbox is not None:
        ax.set_xlim((bbox[0], bbox[1]))
        ax.set_ylim((bbox[2], bbox[3]))

    ax.set_title(title, fontsize=title_size)

    if include_axis:
        ax.set_xlabel("Easting (m)")
        ax.set_ylabel("Northing (m)")
    else:
        ax.axis("off")

    if save_file is not None:
        plt.savefig(save_file, dpi=300, bbox_inches="tight")

    if new_plot:
        plt.show()
    else:
        return ax