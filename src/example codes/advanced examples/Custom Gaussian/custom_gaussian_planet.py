import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import csv

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

rigidity_inputs = {
    "startrigidity": 20,
    "endrigidity": 0,
    "rigiditystep": 0.01,
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
    "intmodel": "Boris-Buneman",
    "gyropercent": 5,
    "minaltitude": 20,
    "maxdistance": 100,
    "maxtime": 0,
    "mintrapdist": 6.6,
    "startaltitude": 20,
    "betaerror": 0.001,
    "totalbetacheck": True,
    "adaptivestep": True,
    "maxsteps": 1000000
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
    "longstep": 5,
    "maxlat": 90,
    "minlat": -90,
    "maxlong": 360,
    "minlong": 0
}


# ============================================================
# RUN OTSO
# ============================================================

def run_otso(glist, hlist):

    print("Running OTSO planet calculation...")

    custom_field_inputs = {
    "g": glist,
    "h": hlist,
    "max_degree": 10,
    "MHDfile": None,
    "MHDcoordsys": None
    }

    magfield_inputs = {
    "internalmag": "Custom Gauss",
    "externalmag": "NONE",
    "boberg": False,
    "bobergtype": "EXTENSION",
    "magnetopause": "Sphere",
    "spheresize": 10,
    "AdaptiveExternalModel": False
}

    planet_results = planet(
        cutoff_comp="Vertical",

        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
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
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs,
        grid_params=grid_inputs
    )

    # Save the OTSO result
    planet_results[0].to_csv("planet.csv", index=False)

    print("OTSO calculation complete.")
    print("Saved results to planet.csv")

    return planet_results[0]


# ============================================================
# PLOT PLANET DATA
# ============================================================

def plot_planet(planet_data):

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

    values = np.arange(0, 5.5, 0.5)

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
        transform=ccrs.PlateCarree()
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
        ticks=np.arange(0, 20, 2)
    )

    color.set_label(
        label="Cut-off Rigidity [GV]",
        size=20
    )

    color.ax.tick_params(
        labelsize=15
    )

    plt.tight_layout()
    plt.title("40950BP Laschamps Cut-offs", size=20)

    plt.savefig(
        "planetplot.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    from pymagglobal import Model, coefficients
    from transform_to_g_and_h import transform_to_g_and_h

    lsmod = Model('LSMOD.2')
    _, _, coeffs = coefficients(1950-40950, lsmod)
    _, _, (glist, hlist) = transform_to_g_and_h(coeffs, lmax_out=13)

    # 1. Run OTSO
    planet_data = run_otso(glist, hlist)

    # 2. Plot the OTSO output
    plot_planet(planet_data)
