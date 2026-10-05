import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter
from scipy.ndimage import uniform_filter


# ============================================================
# SETTINGS
# ============================================================

CSV_PATH = "planet.csv"
OUT_PATH = "planetplot.png"

# Box-filter width (in grid cells) used to smooth out point-to-point noise
# from the MHD-grid field interpolation. 1 disables smoothing; 3 averages
# each point with its 8 immediate neighbours, 5 with its 24 neighbours, etc.
SMOOTH_SIZE = 3


# ============================================================
# LOAD
# ============================================================

def load_grid(csv_path):
    df = pd.read_csv(csv_path)

    lon = np.sort(df["Longitude"].unique())
    lat = np.sort(df["Latitude"].unique())

    Z = (
        df
        .pivot_table(index="Latitude", columns="Longitude", values="Rc [GV]")
        .reindex(index=lat, columns=lon)
        .to_numpy()
    )

    return lon, lat, Z


# ============================================================
# SMOOTH
# ============================================================

def smooth_grid(Z, lon, size):
    """Replace each grid point with the mean of its size x size neighbourhood.
    Latitude is clamped at the poles; longitude wraps around since it's periodic."""

    if size <= 1:
        return Z

    # If both the 0 and 360 longitude columns are present they're the same
    # physical meridian - drop the duplicate 360 edge before wrapping so it
    # isn't double-counted, then restore it as a copy of the smoothed 0 column.
    has_dup_edge = lon[0] == 0 and lon[-1] == 360
    core = Z[:, :-1] if has_dup_edge else Z

    smoothed = uniform_filter(core, size=size, mode=["nearest", "wrap"])

    if has_dup_edge:
        smoothed = np.concatenate([smoothed, smoothed[:, :1]], axis=1)

    # Cut-off rigidity can't be negative; uniform_filter's floating-point
    # noise otherwise nudges near-zero cells (e.g. the polar caps) just
    # below 0, which contourf then leaves unfilled since it's below the
    # lowest plotted level.
    return np.clip(smoothed, 0, None)


# ============================================================
# PLOT
# ============================================================

def plot_planet(lon, lat, Z, out_path=OUT_PATH):

    X, Y = np.meshgrid(lon, lat)

    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig = plt.figure(figsize=(8, 5), dpi=300)

    ax = fig.add_subplot(
        1,
        1,
        1,
        projection=ccrs.Robinson(central_longitude=0)
    )

    ax.set_global()

    # --------------------------------------------------------
    # Gridlines
    # --------------------------------------------------------

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        linewidth=0.8,
        color="black",
        alpha=0.5,
        linestyle="-",
        draw_labels=False
    )

    gl.xlines = True
    gl.ylines = True

    gl.xlocator = mticker.FixedLocator(
        np.arange(-180, 181, 60)
    )

    gl.ylocator = mticker.FixedLocator(
        np.arange(-90, 91, 30)
    )

    def degree_formatter(x, pos):
        return f"{int(x)}°"

    gl.xformatter = FuncFormatter(degree_formatter)
    gl.yformatter = FuncFormatter(degree_formatter)

    gl.xlabel_style = {
        "size": 15
    }

    gl.ylabel_style = {
        "size": 15
    }

    # --------------------------------------------------------
    # Colour scale
    # --------------------------------------------------------

    values = np.arange(0, 20.5, 0.5)

    cmap = plt.get_cmap("viridis")

    norm = BoundaryNorm(
        boundaries=values,
        ncolors=cmap.N,
        clip=True
    )

    sm = ScalarMappable(
        cmap=cmap,
        norm=norm
    )

    sm.set_array([])

    # --------------------------------------------------------
    # Contour plot
    # --------------------------------------------------------

    ax.contourf(
        X,
        Y,
        Z,
        levels=values,
        cmap=cmap,
        norm=norm,
        extend="both",
        transform=ccrs.PlateCarree()
    )

    # --------------------------------------------------------
    # Coastlines
    # --------------------------------------------------------

    ax.coastlines(
        zorder=1
    )

    # --------------------------------------------------------
    # Colour bar
    # --------------------------------------------------------

    color = fig.colorbar(
        sm,
        orientation="horizontal",
        ax=ax,
        fraction=0.07,
        pad=0.08,
        ticks=np.arange(0, 22, 2)
    )

    color.set_label(
        label="Cut-off Rigidity [GV]",
        size=20
    )

    color.ax.tick_params(
        labelsize=15
    )

    # --------------------------------------------------------
    # Final formatting
    # --------------------------------------------------------

    plt.tight_layout()

    plt.savefig(
        out_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    lon, lat, Z = load_grid(CSV_PATH)

    print("Maximum cut-off rigidity:", np.nanmax(Z))
    print("Minimum cut-off rigidity:", np.nanmin(Z))

    Z = smooth_grid(Z, lon, SMOOTH_SIZE)

    plot_planet(lon, lat, Z)
