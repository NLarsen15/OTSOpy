from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from OTSO import skymap

date_inputs           = {"year": 2020, "month": 1, "day": 1, "hour":1,
                        "minute":0, "second":0}

magfield_inputs       = {"internalmag": "IGRF", "externalmag": "TSY89c",
                        "boberg": False, "bobergtype": "EXTENSION",
                        "magnetopause": "Sphere", "spheresize": 10,
                        "AdaptiveExternalModel": False}

rigidity_inputs       = {"startrigidity": 100, "endrigidity": 0,
                        "rigiditystep": 0.01, "rigidityscan": "ON"}

solar_wind_inputs     = {"vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 1, "bz": 1,
                         "by_avg": 0, "bz_avg": 0, "density": 1, "pdyn": 0}

geomagnetic_inputs    = {"Dst": 0, "kp": 2, "n_index": 0, "b_index": 0,
                         "sym_h_corrected": 0}

tsyganenko_inputs     = {"G1": 0, "G2": 0, "G3": 0, "W1": 0, "W2": 0, "W3": 0,
                         "W4": 0, "W5": 0, "W6": 0}

integration_inputs    = {"intmodel": "Boris-Buneman", "gyropercent": 15, "minaltitude": 20,
                         "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                         "startaltitude": 20, "betaerror": 0.001,
                         "totalbetacheck": False, "adaptivestep": False,
                         "maxsteps": 100000}

particle_inputs       = {"Anum": 1, "anti": "YES"}

computation_inputs    = {"corenum": 5, "threadnum": 1, "Verbose": True}

custom_field_inputs   = {"g": None, "h": None, "max_degree": 13, "MHDfile": None,
                         "MHDcoordsys": None}

coord_inputs          = {"inputcoord": "GDZ"}
# Note: coordsystem will also control the coordinate system of asymptotic viewing direction results

data_retrieval_inputs = {"serverdata": "OFF", "livedata": "OFF"}

skymap_inputs         = {"zenithstep": 15, "azimuthstep": 45, "maxzenith": 75,
                         "minzenith":0, "maxazimuth": 360, "minazimuth": 0}

def plot_skymap(df, value="Rc [GV]", title=None):

    df = df.copy()

    zenith_zero = df[df["Zenith"] == 0]

    if len(zenith_zero) == 1:

        base = zenith_zero.iloc[0]

        fill_rows = []

        for az in np.arange(0, 360, 45):
            row = base.copy()
            row["Azimuth"] = az
            fill_rows.append(row)

        df = pd.concat(
            [
                df[df["Zenith"] != 0],
                pd.DataFrame(fill_rows)
            ],
            ignore_index=True
        )

    df.loc[df["Azimuth"] == 360, "Azimuth"] = 0

    theta = np.deg2rad(df["Azimuth"].values)
    radius = df["Zenith"].values
    values = df[value].values

    theta = np.concatenate([theta, theta + 2*np.pi, theta - 2*np.pi])
    radius = np.concatenate([radius, radius, radius])
    values = np.concatenate([values, values, values])

    theta_grid, radius_grid = np.meshgrid(
        np.linspace(0, 2*np.pi, 361),
        np.linspace(0, radius.max(), 200)
    )

    Z = griddata(
        (theta, radius),
        values,
        (theta_grid, radius_grid),
        method="linear"
    )

    if np.isnan(Z).any():

        Z2 = griddata(
            (theta, radius),
            values,
            (theta_grid, radius_grid),
            method="nearest"
        )

        Z[np.isnan(Z)] = Z2[np.isnan(Z)]

    fig = plt.figure(figsize=(9, 9))
    ax = plt.subplot(111, projection="polar")

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    contour = ax.contourf(
        theta_grid,
        radius_grid,
        Z,
        levels=30,
        cmap="viridis"
    )

    ax.grid(
        True,
        linewidth=1.2,
        alpha=0.7
    )

    ax.set_ylim(0, radius.max())

    ax.set_yticks(
        [0, 15, 30, 45, 60, 75]
    )

    ax.set_yticklabels(
        ["0°", "15°", "30°", "45°", "60°", "75°"],
        fontsize=13,
        fontweight="bold"
    )

    ax.tick_params(
        axis="y",
        pad=10,
        width=1.5,
        length=6
    )

    ax.tick_params(
        axis="x",
        labelsize=13,
        width=1.5,
        length=6
    )

    for label in ax.get_xticklabels():
        label.set_fontweight("bold")

    cbar = plt.colorbar(
        contour,
        ax=ax,
        pad=0.12,
        shrink=0.8
    )

    cbar.set_label(
        value,
        fontsize=14,
        fontweight="bold"
    )

    cbar.ax.tick_params(
        labelsize=12
    )

    plt.tight_layout()
    plt.savefig("skymap.png", dpi=400)
    plt.show()


if __name__ == "__main__":

    stations = ["ROME"]

    dfs = {}

    skymap_results = skymap(
            Stations=stations,
            skymap_params=skymap_inputs,
            datetime_params=date_inputs,
            magfield_params=magfield_inputs,
            rigidity_params=rigidity_inputs,
            solar_wind_params=solar_wind_inputs,
            geomagnetic_params=geomagnetic_inputs,
            tsyganenko_params=tsyganenko_inputs,
            integration_params=integration_inputs,
            particle_params=particle_inputs,
            computation_params=computation_inputs,
            coordinate_params=coord_inputs,
            custom_field_params=custom_field_inputs,
            data_retrieval_params=data_retrieval_inputs
        )

    # skymap_results[0] is the dictionary of DataFrames
    result_dict = skymap_results[0]

    for station in stations:

            dfs[station] = result_dict[station]

            csv_file = Path(f"{station}_skymap.csv")

            dfs[station].to_csv(csv_file, index=False)

            print(f"Saved {csv_file}")

    # --------------------------------------------------
    # Plot
    # --------------------------------------------------

    plot_skymap(
        dfs["ROME"],
        value="Rc [GV]",
        title="ROME Cutoff Rigidity"
    )