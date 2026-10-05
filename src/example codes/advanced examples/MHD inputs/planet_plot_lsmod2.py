import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

from matplotlib.colors import BoundaryNorm, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter

from pymagglobal import Model, coefficients
from transform_to_g_and_h import transform_to_g_and_h

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
    "gyropercent": 1,
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

# The LSMOD.2 g/h coefficients used both to build LSMOD2_grid_uniform_01.csv
# (gen_lsmod2_uniform.py) and, here, as the ungridded "baseline" reference -
# same epoch/lmax_out/max_degree as the generator, so the baseline is an
# exact evaluation of the same field the grid approximates.
_epoch = 40950  # BP
_lsmod = Model('LSMOD.2')
_, _, _coeffs = coefficients(1950 - _epoch, _lsmod)
_, _, (_glist, _hlist) = transform_to_g_and_h(_coeffs, lmax_out=13)

custom_field_inputs = {"g": None, "h": None, "max_degree": 13, "MHDcoordsys": "GEO",
                        "MHDfile": "LSMOD2_grid_uniform_01.csv"}

# Same field, but the geometrically-stretched grid from
# gen_lsmod2_stretched.py - 0.001 Re spacing anchored right at Earth's
# surface (x=+-1 Re), coarsening both inward (unused interior) and outward
# to +/-6.6 Re - instead of the fixed uniform grid. MHDgridtype must be
# "stretched" since the spacing isn't constant. Generated separately; run
# gen_lsmod2_stretched.py first if this file doesn't exist yet.
custom_field_inputs_stretched = {"g": None, "h": None, "max_degree": 13, "MHDcoordsys": "GEO",
                                  "MHDfile": "LSMOD2_grid_stretched_r1_0001_ext66.csv",
                                  "MHDgridtype": "stretched"}

# Runs to compare: the LSMOD.2 paleomagnetic field (40,950 BP) interpolated
# four ways - "trilinear" (the old default), "tricubic" (smoother, can
# overshoot near sharp gradients), "monotonic" (same tricubic fit, but
# limited so it never overshoots), and "divfree" (interpolates the vector
# potential and returns B = curl(A), divergence-free by construction; see
# CustomFieldParams docs for details) - on both the uniform grid and the
# surface-refined stretched grid, plus a "normal" baseline that swaps the
# MHD grid for a direct (ungridded) evaluation of the same LSMOD.2 g/h
# coefficients via internalmag="Custom Gauss", keeping every other magfield
# parameter identical, so the only thing that differs is whether the field
# goes through the grid/interpolation path at all. "divfree" requires a
# uniform grid (the FFT-based vector potential solve assumes constant grid
# spacing), so there's no stretched-grid counterpart for it. Stretched-grid
# runs get their own labels (not just an overwrite of the uniform-grid ones)
# so both sets of results - and the skip-if-CSV-exists cache in the main
# loop below - stay independent.
RUNS = [
    {
        "label": "lsmod2_trilinear",
        "title": "LSMOD.2 interpolation: trilinear (uniform grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs, "MHDinterpolation": "trilinear"},
    },
    {
        "label": "lsmod2_tricubic",
        "title": "LSMOD.2 interpolation: tricubic (uniform grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs, "MHDinterpolation": "tricubic"},
    },
    {
        "label": "lsmod2_monotonic",
        "title": "LSMOD.2 interpolation: monotonic (uniform grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs, "MHDinterpolation": "monotonic"},
    },
    {
        "label": "lsmod2_divfree",
        "title": "LSMOD.2 interpolation: divfree (uniform grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs, "MHDinterpolation": "divfree"},
    },
    {
        "label": "lsmod2_trilinear_stretched",
        "title": "LSMOD.2 interpolation: trilinear (stretched grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs_stretched, "MHDinterpolation": "trilinear"},
    },
    {
        "label": "lsmod2_tricubic_stretched",
        "title": "LSMOD.2 interpolation: tricubic (stretched grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs_stretched, "MHDinterpolation": "tricubic"},
    },
    {
        "label": "lsmod2_monotonic_stretched",
        "title": "LSMOD.2 interpolation: monotonic (stretched grid)",
        "magfield_params": magfield_inputs,
        "custom_field_params": {**custom_field_inputs_stretched, "MHDinterpolation": "monotonic"},
    },
    {
        "label": "lsmod2_baseline",
        "title": "Baseline: LSMOD.2 Custom Gauss (no MHD grid)",
        "magfield_params": {**magfield_inputs, "internalmag": "Custom Gauss", "externalmag": "NONE"},
        "custom_field_params": {"g": _glist, "h": _hlist, "max_degree": 10,
                                 "MHDfile": None, "MHDcoordsys": None},
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

    values = np.arange(0, 5.25, 0.25)

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
        # Project the coordinate grid first, then let matplotlib contour the
        # already-projected data natively - avoids cartopy's polygon-based
        # reprojection path, which can crash on a Robinson projection with
        # "'GeometryCollection' object is not subscriptable" when a filled
        # band touches the map edge or the data contains NaNs (both come up
        # here: NaN cutoffs are a real possibility for a paleomagnetic field
        # like LSMOD.2 that may be very weak or non-dipolar at some epochs).
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
        ticks=np.arange(0, 5.5, 0.5)
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


# ============================================================
# DIFFERENCE PLOTS (interpolation method vs. baseline)
# ============================================================

def to_grid(planet_data):
    """Pivot a planet() result's (Latitude, Longitude, Rc) rows onto a
    regular lat/lon grid, same layout plot_planet() builds internally."""
    lat = np.sort(np.asarray(planet_data["Latitude"], dtype=float))
    lat = np.unique(lat)
    lon = np.sort(np.unique(np.asarray(planet_data["Longitude"], dtype=float)))
    Z = (
        planet_data
        .assign(Longitude=planet_data["Longitude"].astype(float),
                Latitude=planet_data["Latitude"].astype(float))
        .pivot_table(index="Latitude", columns="Longitude", values="Rc [GV]")
        .reindex(index=lat, columns=lon)
        .to_numpy()
    )
    return lat, lon, Z


def plot_diff(lat, lon, diffZ, label, title, vmax):

    X, Y = np.meshgrid(lon, lat)

    fig = plt.figure(figsize=(8, 5), dpi=300)

    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
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
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
    gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 30))

    def degree_formatter(x, pos):
        return f"{int(x)}°"

    gl.xformatter = FuncFormatter(degree_formatter)
    gl.yformatter = FuncFormatter(degree_formatter)
    gl.xlabel_style = {"size": 15}
    gl.ylabel_style = {"size": 15}

    cmap = plt.get_cmap("RdBu_r")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    levels = np.linspace(-vmax, vmax, 21)

    ax.contourf(
        X, Y, diffZ,
        levels=levels,
        cmap=cmap,
        norm=norm,
        extend="both",
        transform=ccrs.PlateCarree(),
        transform_first=True  # see plot_planet() for why
    )

    ax.coastlines(zorder=1)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    color = fig.colorbar(
        sm,
        orientation="horizontal",
        ax=ax,
        fraction=0.07,
        pad=0.08,
    )
    color.set_label(label="Δ Cut-off Rigidity [GV] (method − baseline)", size=18)
    color.ax.tick_params(labelsize=15)

    plt.tight_layout()

    out_png = f"planetplot_diff_{label}.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved diff plot to {out_png}")

    plt.close(fig)


if __name__ == "__main__":

    results = {}

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

        results[run["label"]] = planet_data

    # 3. Difference plots: each interpolation method vs. the ungridded
    # baseline. All three share one symmetric color scale (the largest
    # |diff| across all of them) so the maps are directly comparable to
    # each other, not just internally self-consistent.
    baseline_runs = [r for r in RUNS if r["label"].endswith("_baseline")]
    method_runs = [r for r in RUNS if not r["label"].endswith("_baseline")]

    if baseline_runs:
        baseline_label = baseline_runs[0]["label"]
        base_lat, base_lon, base_Z = to_grid(results[baseline_label])

        diffs = {}
        for run in method_runs:
            lat, lon, Z = to_grid(results[run["label"]])
            if not (np.array_equal(lat, base_lat) and np.array_equal(lon, base_lon)):
                print(f"Skipping diff for {run['label']}: lat/lon grid doesn't match baseline")
                continue
            diffs[run["label"]] = (lat, lon, Z - base_Z)

        if diffs:
            vmax = max(np.nanmax(np.abs(d)) for _, _, d in diffs.values())
            print(f"\nDifference plots vs. {baseline_label}: shared scale +/-{vmax:.3f} GV")
            for run in method_runs:
                if run["label"] not in diffs:
                    continue
                lat, lon, diffZ = diffs[run["label"]]
                title = f"{run['title']} minus baseline"
                plot_diff(lat, lon, diffZ, run["label"], title, vmax)
    else:
        print("\nNo run labeled '*_baseline' found - skipping difference plots.")
