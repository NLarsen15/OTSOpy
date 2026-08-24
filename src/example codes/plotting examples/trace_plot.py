from OTSO import trace
import matplotlib.pyplot as plt
import matplotlib.patches as patches


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
    "externalmag": "NONE",
    "boberg": False,
    "bobergtype": "EXTENSION",
    "magnetopause": "NONE"
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

grid_inputs = {
    "latstep": -5,
    "longstep": 180,
    "maxlat": 85,
    "minlat": -85,
    "maxlong": 360,
    "minlong": 0
}

integration_inputs = {
    "minaltitude": 20,
    "maxdistance": 100,
    "maxtime": 0,
    "startaltitude": 20,
    "maxsteps": 100000
}

computation_inputs = {
    "corenum": 4,
    "Verbose": True
}

custom_field_inputs = {
    "g": None,
    "h": None,
    "max_degree": 13,
    "MHDfile": None,
    "MHDcoordsys": None
}

data_retrieval_inputs = {
    "serverdata": "OFF",
    "livedata": "OFF"
}


def plot_xz(trace_results):

    # trace_results[0] contains the dictionary of traces
    traces = trace_results[0]

    # Wider figure
    fig, ax = plt.subplots(figsize=(14, 8))

    # Plot all field lines
    for key, data in traces.items():

        trace_df = data["Trace"]

    # Keep only points within the requested Z range
        mask = (
            (trace_df["Y_GEO [Re]"] > -2) &
            (trace_df["Y_GEO [Re]"] < 2)
        )

        # Keep original indexing and insert NaNs outside the Y range
        x = trace_df["X_GEO [Re]"].where(mask)
        z = trace_df["Z_GEO [Re]"].where(mask)

        ax.plot(
            x,
            z,
            linewidth=1.5,
            color="black"
        )

    # Earth
    earth = patches.Circle(
        (0, 0),
        1,
        facecolor="blue",
        edgecolor="black",
        linewidth=2,
        alpha=0.8
    )

    ax.add_patch(earth)

    # Labels
    ax.set_xlabel(
        "X (GEO) [Re]",
        fontsize=18
    )

    ax.set_ylabel(
        "Z (GEO) [Re]",
        fontsize=18
    )

    # Bigger tick labels
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=16
    )

    # Keep Earth circular
    ax.set_aspect("equal", adjustable="box")

    # ---------------------------------------------------------
    # VIEWING WINDOW
    #
    # IMPORTANT:
    # These limits DO NOT delete/cut the field-line data.
    # They only determine what part of the plot is visible.
    # ---------------------------------------------------------

    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)

    # Optional grid
    ax.grid(
        True,
        alpha=0.3
    )

    plt.tight_layout()

    plt.savefig(
        "Trace_Plot.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


if __name__ == '__main__':

    coordinatesystem = "GEO"

    trace_results = trace(
        coordsys=coordinatesystem,
        datetime_params=date_inputs,
        computation_params=computation_inputs,
        magfield_params=magfield_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        grid_params=grid_inputs
    )

    print(trace_results[0])

    plot_xz(trace_results)