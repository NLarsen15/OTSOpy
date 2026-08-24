import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D

from OTSO import trajectory


# ============================================================
# Plotting utilities
# ============================================================

def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)

    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)

    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


# ============================================================
# OTSO input parameters
# ============================================================

date_inputs = {
    "year": 2020,
    "month": 1,
    "day": 1,
    "hour": 1,
    "minute": 0,
    "second": 0
}

magfield_inputs = {
    "internalmag": "IGRF",
    "externalmag": "TSY89c",
    "boberg": False,
    "bobergtype": "EXTENSION",
    "magnetopause": "Sphere",
    "spheresize": 10,
    "AdaptiveExternalModel": False
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
    "adaptivestep": False,
    "maxsteps": 1000000
}

particle_inputs = {
    "rigidity": 0.4,
    "Anum": 1,
    "anti": "YES",
    "zenith": 0,
    "azimuth": 0
}

computation_inputs = {
    "corenum": 1,
    "Verbose": True
}

custom_field_inputs = {
    "g": None,
    "h": None,
    "max_degree": 13,
    "MHDfile": None,
    "MHDcoordsys": None
}

coord_inputs = {
    "coordsystem": "GEO",
    "inputcoord": "GDZ"
}

data_retrieval_inputs = {
    "serverdata": "OFF",
    "livedata": "OFF"
}


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    # Neutron monitor stations
    stations_list = [
        "OULU",
    ]

    # --------------------------------------------------------
    # Run OTSO
    # --------------------------------------------------------

    trajectory_results = trajectory(
        Stations=stations_list,

        particle_params=particle_inputs,
        computation_params=computation_inputs,

        datetime_params=date_inputs,

        magfield_params=magfield_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,

        integration_params=integration_inputs,

        custom_field_params=custom_field_inputs,
        coordinate_params=coord_inputs,

        data_retrieval_params=data_retrieval_inputs,
    )

    # --------------------------------------------------------
    # Inspect OTSO output
    # --------------------------------------------------------

    print("First OTSO result:")
    print(trajectory_results[0])

    print("\nInput information:")
    print(trajectory_results[-1])


    # ========================================================
    # Load Earth texture
    # ========================================================

    texture = mpimg.imread("Blue_Marble.jpg")

    downsample_factor = 40
    texture_ds = texture[::downsample_factor, ::downsample_factor]

    rows, cols = texture_ds.shape[:2]

    u = np.linspace(0, 2 * np.pi, cols) + np.pi
    v = np.linspace(0, np.pi, rows)

    u, v = np.meshgrid(u, v)

    earth_radius = 1.0

    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)


    # ========================================================
    # Create 3D plot
    # ========================================================

    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(
        111,
        projection="3d"
    )


    # ========================================================
    # Plot Earth
    # ========================================================

    if texture_ds.max() > 1:
        texture_colors = texture_ds / 255.0
    else:
        texture_colors = texture_ds

    ax.plot_surface(
        x,
        y,
        z,
        rstride=1,
        cstride=1,
        facecolors=texture_colors,
        linewidth=0,
        antialiased=False
    )


    # ========================================================
    # Plot OTSO trajectories
    # ========================================================

    trajectory_list = trajectory_results[0]

    for result in trajectory_list:

        station = result["station"]

        # DataFrame containing the trajectory
        df = result["trajectory"]

        # Position in GEO coordinates [Re]
        x_traj = df["GEO_X [Re]"].to_numpy()
        y_traj = df["GEO_Y [Re]"].to_numpy()
        z_traj = df["GEO_Z [Re]"].to_numpy()

        ax.plot(
            x_traj,
            y_traj,
            z_traj,
            linewidth=1.5,
            label=f"{station} (R={result['rigidity']:.2f} GV)"
        )


    # ========================================================
    # Plot formatting
    # ========================================================

    ax.set_xlabel("X [Re]")
    ax.set_ylabel("Y [Re]")
    ax.set_zlabel("Z [Re]")

    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([-5, 5])

    ax.set_box_aspect([1, 1, 1])

    ax.set_title(
        "OTSO Particle Trajectory"
    )

    ax.legend()

    plt.tight_layout()
    plt.show()