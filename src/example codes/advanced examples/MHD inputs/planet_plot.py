import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter

from OTSO import planet


# ============================================================
# OTSO INPUT PARAMETERS
# ============================================================

date_inputs = {
    "year": 2020,
    "month": 1,
    "day": 1,
    "hour": 1,
    "minute": 0,
    "second": 0
}

magfield_inputs       = {"internalmag": "NONE", "externalmag": "MHD",
                        "boberg": False, "bobergtype": "EXTENSION",
                        "magnetopause": "Sphere", "spheresize": 6.6,
                        "AdaptiveExternalModel": False}

rigidity_inputs = {
    "startrigidity": 20,
    "endrigidity": 0,
    "rigiditystep": 0.1,
    "rigidityscan": "ON"
}

asymptotic_inputs = {
    "asymptotic": "NO",
    "asymlevels": [
        0.1, 0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9,
        10, 15, 20, 30, 50, 70, 100, 300, 500, 700, 1000
    ],
    "unit": "GeV"
}

transmission_inputs = {
    "transmission": False,
    "transmissionRstep": 0.001,
    "transmissionsamples": 20
}

solar_wind_inputs = {
    "vx": -500,
    "vy": 0,
    "vz": 0,
    "bx": 0,
    "by": 1,
    "bz": 1,
    "by_avg": 0,
    "bz_avg": 0,
    "density": 1,
    "pdyn": 0
}

geomagnetic_inputs = {
    "Dst": 0,
    "kp": 2,
    "n_index": 0,
    "b_index": 0,
    "sym_h_corrected": 0
}

tsyganenko_inputs = {
    "G1": 0,
    "G2": 0,
    "G3": 0,
    "W1": 0,
    "W2": 0,
    "W3": 0,
    "W4": 0,
    "W5": 0,
    "W6": 0
}

integration_inputs = {
    "intmodel": "6RK",
    "gyropercent": 15,
    "minaltitude": 20,
    "maxdistance": 100,
    "maxtime": 0,
    "mintrapdist": 6.6,
    "startaltitude": 20,
    "betaerror": 0.01,
    "totalbetacheck": False,
    "adaptivestep": True,
    "maxsteps": 100000
}

particle_inputs = {
    "Anum": 1,
    "anti": "YES",
    "zenith": 0,
    "azimuth": 0
}

computation_inputs = {
    "corenum": 7,
    "threadnum": 1,
    "Verbose": True,
    "delim": ";"
}

custom_field_inputs   = {"g": None, "h": None, "max_degree": 13, "MHDcoordsys": "GEO",
                          "MHDfile": "IGRF_grid_uniform_01.csv"}

# Same grid, but the geometrically-stretched one from gen_igrf_stretched.py -
# 0.001 Re spacing anchored right at Earth's surface (x=+-1 Re), coarsening
# both inward (unused interior) and outward to +/-6.6 Re (matching the
# magnetopause spheresize/mintrapdist used below) - instead of the fixed
# 0.15 Re uniform grid. MHDgridtype must be "stretched" since the spacing
# isn't constant. Generated separately; run gen_igrf_stretched.py first if
# this file doesn't exist yet.
custom_field_inputs_stretched = {"g": None, "h": None, "max_degree": 13, "MHDcoordsys": "GEO",
                                  "MHDfile": "IGRF_grid_stretched_r1_0001_ext66.csv",
                                  "MHDgridtype": "stretched"}

# Runs to compare: the MHD field interpolated three ways - "trilinear" (the
# old default), "tricubic" (smoother, can overshoot near sharp gradients),
# and "monotonic" (same tricubic fit, but limited so it never overshoots;
# see CustomFieldParams docs for details) - on both the uniform grid and the
# near-Earth-refined stretched grid, plus a "normal" baseline that swaps the
# MHD external field for plain IGRF (internalmag="IGRF", externalmag="NONE"),
# keeping every other magfield parameter (magnetopause, spheresize, etc.)
# identical to the MHD runs, so the only thing that differs is the field
# source itself. Stretched-grid runs get their own labels (not just an
# overwrite of the uniform-grid ones) so both sets of results - and the
# skip-if-CSV-exists cache in the main loop below - stay independent.
RUNS = [
    {
        "label": "trilinear",
        "title": "MHD interpolation: trilinear (uniform grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs, "MHDinterpolation": "trilinear"},
    },
    {
        "label": "tricubic",
        "title": "MHD interpolation: tricubic (uniform grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs, "MHDinterpolation": "tricubic"},
    },
    {
        "label": "monotonic",
        "title": "MHD interpolation: monotonic (uniform grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs, "MHDinterpolation": "monotonic"},
    },
    {
        "label": "trilinear_stretched",
        "title": "MHD interpolation: trilinear (stretched grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs_stretched, "MHDinterpolation": "trilinear"},
    },
    {
        "label": "tricubic_stretched",
        "title": "MHD interpolation: tricubic (stretched grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs_stretched, "MHDinterpolation": "tricubic"},
    },
    {
        "label": "monotonic_stretched",
        "title": "MHD interpolation: monotonic (stretched grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs_stretched, "MHDinterpolation": "monotonic"},
    },
    {
        "label": "igrf_baseline",
        "title": "Baseline: IGRF only (no MHD grid)",
        "magfield_params": {**magfield_inputs, "internalmag": "IGRF", "externalmag": "NONE"},
        "custom_field_params": {},
    },
]

coord_inputs = {
    "coordsystem": "GEO",
    "inputcoord": "GDZ"
}

data_retrieval_inputs = {
    "serverdata": "OFF",
    "livedata": "OFF"
}

grid_inputs = {
    "latstep": -5,
    "longstep": 15,
    "maxlat": 90,
    "minlat": -90,
    "maxlong": 360,
    "minlong": 0
}


# ============================================================
# RUN OTSO
# ============================================================

def csv_path_for(label):
    return f"planet_{label}.csv"


def run_otso(run):

    print(f"Running OTSO planet calculation ({run['label']})...")

    planet_results = planet(
        cutoff_comp="Vertical",

        datetime_params=date_inputs,
        magfield_params=run["magfield_params"],
        rigidity_params=rigidity_inputs,
        asymptotic_params=asymptotic_inputs,
        transmission_params=transmission_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params=integration_inputs,
        particle_params=particle_inputs,
        computation_params=computation_inputs,
        coordinate_params=coord_inputs,
        custom_field_params=run["custom_field_params"],
        data_retrieval_params=data_retrieval_inputs,
        grid_params=grid_inputs
    )

    # Save the OTSO result, one CSV per run so they don't overwrite each other.
    out_csv = csv_path_for(run["label"])
    planet_results[0].to_csv(out_csv, index=False)

    print("OTSO calculation complete.")
    print(f"Saved results to {out_csv}")

    return planet_results[0]


# ============================================================
# PLOT PLANET DATA
# ============================================================

def plot_planet(planet_data, label, title):


    PlanetLat = np.asarray(planet_data["Latitude"], dtype=float)
    PlanetLong = np.asarray(planet_data["Longitude"], dtype=float)
    PlanetDose = np.asarray(planet_data["Rc [GV]"], dtype=float)

    print("Maximum cut-off rigidity:", np.nanmax(PlanetDose))
    print("Minimum cut-off rigidity:", np.nanmin(PlanetDose))


    PlanetDf = pd.DataFrame({
        "x": PlanetLong,
        "y": PlanetLat,
        "z": PlanetDose
    })


    X_unique = np.sort(PlanetDf["x"].unique())
    Y_unique = np.sort(PlanetDf["y"].unique())

    X, Y = np.meshgrid(X_unique, Y_unique)

    Z = (
        PlanetDf
        .pivot_table(
            index="y",
            columns="x",
            values="z"
        )
        .reindex(
            index=Y_unique,
            columns=X_unique
        )
        .values
    )

    fig = plt.figure(figsize=(8, 5), dpi=300)

    ax = fig.add_subplot(
        1,
        1,
        1,
        projection=ccrs.Robinson(central_longitude=0)
    )

    ax.set_global()

    ax.set_title(title, fontsize=16)


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

    gl.ylocator = mticker.FixedLocator(
        np.arange(-90, 91, 30)
    )

    gl.xlocator = mticker.FixedLocator(
        np.arange(-180, 181, 60)
    )

    gl.xlabel_style = {
        "size": 15
    }

    gl.ylabel_style = {
        "size": 15
    }

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

    plot = ax.contourf(
        X,
        Y,
        Z,
        levels=values,
        cmap=cmap,
        norm=norm,
        transform=ccrs.PlateCarree(),
        # Project the grid first, then contour natively in projected space -
        # avoids a cartopy/shapely bug where reprojecting filled contour
        # polygons on a Robinson projection can raise "'GeometryCollection'
        # object is not subscriptable" (seen in planet_plot_lsmod2.py, e.g.
        # from a NaN cutoff or a filled band touching the map edge).
        transform_first=True
    )

    ax.coastlines(
        zorder=1
    )

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

    plt.tight_layout()

    out_png = f"planetplot_{label}.png"

    plt.savefig(
        out_png,
        dpi=300,
        bbox_inches="tight"
    )

    print(f"Saved plot to {out_png}")

    # Batch-generating one plot per run - close each figure instead of
    # blocking on plt.show() so the loop doesn't stall waiting for a window
    # to be closed by hand between runs.
    plt.close(fig)

if __name__ == "__main__":

    for run in RUNS:
        print(f"\n=== {run['title']} ===")

        # 1. Run OTSO, unless this run's CSV is already there from a previous
        # pass - the planet cutoff scan is the expensive part, so a rerun
        # (e.g. after tweaking the plot) shouldn't redo it for runs that
        # already finished. Delete the CSV (or the whole file) to force a
        # recompute for a given run.
        out_csv = csv_path_for(run["label"])
        if os.path.exists(out_csv):
            print(f"Found existing {out_csv}, skipping computation and loading cached results.")
            planet_data = pd.read_csv(out_csv)
        else:
            planet_data = run_otso(run)

        # 2. Plot the OTSO output
        plot_planet(planet_data, run["label"], run["title"])
