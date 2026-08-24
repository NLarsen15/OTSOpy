import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

from OTSO import cone


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

magfield_inputs = {
    "internalmag": "IGRF",
    "externalmag": "TSY89c",
    "boberg": False,
    "bobergtype": "EXTENSION",
    "magnetopause": "Sphere",
    "spheresize": 10,
    "AdaptiveExternalModel": False
}

rigidity_inputs = {
    "startrigidity": 20,
    "endrigidity": 0,
    "rigiditystep": 0.001
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
    "gyropercent": 15,
    "minaltitude": 20,
    "maxdistance": 100,
    "maxtime": 0,
    "mintrapdist": 6.6,
    "startaltitude": 20,
    "betaerror": 0.001,
    "totalbetacheck": True,
    "adaptivestep": False,
    "maxsteps": 100000
}

particle_inputs = {
    "Anum": 1,
    "anti": "YES",
    "zenith": 0,
    "azimuth": 0
}

computation_inputs = {
    "corenum": 1,
    "threadnum": 7,
    "Verbose": True,
    "delim": ";"
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
# EXTRACT FILTER FROM OTSO RESULT
# ============================================================

def get_filter(value):

    """
    OTSO returns strings such as:

        -3;-1.6206;241.3022
        -1;56.0704;30.2979
         1;21.3401;273.4899

    The first number is the OTSO filter.

        1  = allowed
       -1  = forbidden
       -3  = forbidden
    """

    if pd.isna(value):
        return None

    text = str(value).strip()

    if text == "":
        return None

    first_field = text.split(";")[0].strip()

    try:
        return int(float(first_field))
    except (ValueError, TypeError):
        return None


# ============================================================
# GET CUTOFFS
# ============================================================

def get_cutoffs(cutoff_df, station):

    station = str(station).strip().upper()

    station_columns = {
        str(column).strip().upper(): column
        for column in cutoff_df.columns
    }

    if station not in station_columns:

        raise ValueError(
            f"Station '{station}' not found.\n"
            f"Available stations: {list(cutoff_df.columns)}"
        )

    actual_station_column = station_columns[station]

    index_map = {
        str(index).strip().upper(): index
        for index in cutoff_df.index
    }

    Rc = float(
        cutoff_df.loc[
            index_map["RC"],
            actual_station_column
        ]
    )

    Ru = float(
        cutoff_df.loc[
            index_map["RU"],
            actual_station_column
        ]
    )

    Rl = float(
        cutoff_df.loc[
            index_map["RL"],
            actual_station_column
        ]
    )

    return Rc, Ru, Rl


# ============================================================
# GET RIGIDITY
# ============================================================

def get_rigidity_values(cone_df):

    """
    IMPORTANT:

    OTSO gives rigidity in the column:

        R [GV]

    Do NOT use the DataFrame index.
    """

    # Exact OTSO column
    if "R [GV]" in cone_df.columns:

        return pd.to_numeric(
            cone_df["R [GV]"],
            errors="coerce"
        )

    # Fallback in case OTSO changes capitalization
    for column in cone_df.columns:

        name = str(column).strip().upper()

        if name == "R [GV]":

            return pd.to_numeric(
                cone_df[column],
                errors="coerce"
            )

    raise ValueError(
        "Could not find the OTSO rigidity column 'R [GV]'.\n"
        f"Available columns: {list(cone_df.columns)}"
    )


# ============================================================
# PREPARE PLOT DATA
# ============================================================

def prepare_plot_data(cone_df, station):

    station = str(station).strip().upper()

    # --------------------------------------------------------
    # Find station column
    # --------------------------------------------------------

    station_columns = {
        str(column).strip().upper(): column
        for column in cone_df.columns
    }

    if station not in station_columns:

        raise ValueError(
            f"Station '{station}' not found.\n"
            f"Available columns: {list(cone_df.columns)}"
        )

    actual_station_column = station_columns[station]

    # --------------------------------------------------------
    # GET ACTUAL RIGIDITY
    # --------------------------------------------------------

    rigidity = get_rigidity_values(
        cone_df
    ).reset_index(drop=True)

    # --------------------------------------------------------
    # GET OTSO RESULT
    # --------------------------------------------------------

    result = (
        cone_df[actual_station_column]
        .reset_index(drop=True)
    )

    # --------------------------------------------------------
    # CREATE DATAFRAME
    # --------------------------------------------------------

    df = pd.DataFrame({
        "R": rigidity,
        "Raw": result
    })

    # --------------------------------------------------------
    # EXTRACT FILTER
    # --------------------------------------------------------

    df["Filter"] = (
        df["Raw"]
        .apply(get_filter)
    )

    # --------------------------------------------------------
    # REMOVE INVALID VALUES
    # --------------------------------------------------------

    df = df.dropna(
        subset=[
            "R",
            "Filter"
        ]
    ).copy()

    df["R"] = df["R"].astype(float)
    df["Filter"] = df["Filter"].astype(int)

    # --------------------------------------------------------
    # ALLOWED = ONLY FILTER 1
    # --------------------------------------------------------

    df["Allowed"] = (
        df["Filter"] == 1
    )

    return df


# ============================================================
# PLOT
# ============================================================

def plot_cone(
    cone_df,
    cutoff_df,
    station,
    xmin=0.50,
    xmax=1.20
):

    station = str(station).strip().upper()

    # ========================================================
    # PREPARE DATA
    # ========================================================

    df = prepare_plot_data(
        cone_df,
        station
    )

    # ========================================================
    # CUTOFFS
    # ========================================================

    Rc, Ru, Rl = get_cutoffs(
        cutoff_df,
        station
    )

    # ========================================================
    # PRINT DIAGNOSTICS
    # ========================================================

    print()
    print("=" * 70)
    print(f"                         {station}")
    print("=" * 70)

    print()
    print(f"Rl = {Rl:.3f} GV")
    print(f"Rc = {Rc:.3f} GV")
    print(f"Ru = {Ru:.3f} GV")

    # --------------------------------------------------------
    # ALL FILTER COUNTS
    # --------------------------------------------------------

    print()
    print("ALL FILTER VALUES:")

    print(
        df["Filter"]
        .value_counts()
        .sort_index()
    )

    # --------------------------------------------------------
    # PLOT RANGE
    # --------------------------------------------------------

    plot_df = df[
        (df["R"] >= xmin) &
        (df["R"] <= xmax)
    ].copy()

    print()
    print(
        f"FILTER VALUES BETWEEN "
        f"{xmin:.2f} AND {xmax:.2f} GV:"
    )

    print(
        plot_df["Filter"]
        .value_counts()
        .sort_index()
    )

    # ========================================================
    # PRINT EVERY POINT IN PLOT RANGE
    # ========================================================

    print()
    print("OTSO POINTS USED FOR PLOTTING:")

    print(
        plot_df[
            ["R", "Filter", "Raw"]
        ].to_string(
            index=False
        )
    )

    # ========================================================
    # CREATE FIGURE
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(20, 6)
    )

    ax.set_facecolor(
        "white"
    )

    fig.patch.set_facecolor(
        "white"
    )

    # ========================================================
    # RECTANGLE WIDTH
    # ========================================================

    rectangle_width = 0.001

    half_width = (
        rectangle_width / 2
    )

    # ========================================================
    # DRAW ONE RECTANGLE PER OTSO POINT
    # ========================================================

    for _, row in plot_df.iterrows():

        R = float(
            row["R"]
        )

        filter_value = int(
            row["Filter"]
        )

        # ----------------------------------------------------
        # FILTER 1 = WHITE
        # ----------------------------------------------------

        if filter_value == 1:

            facecolor = "white"

        # ----------------------------------------------------
        # -1 / -3 = BLACK
        # ----------------------------------------------------

        else:

            facecolor = "black"

        # ----------------------------------------------------
        # Rectangle centered on R
        # ----------------------------------------------------

        rectangle = Rectangle(
            (
                R - half_width,
                0
            ),
            rectangle_width,
            1,
            facecolor=facecolor,
            edgecolor="none",
            linewidth=0,
            zorder=2
        )

        ax.add_patch(
            rectangle
        )

    # ========================================================
    # AXIS LIMITS
    # ========================================================

    ax.set_xlim(
        xmin,
        xmax
    )

    ax.set_ylim(
        0,
        1
    )

    # ========================================================
    # REMOVE Y AXIS
    # ========================================================

    ax.set_yticks([])

    # ========================================================
    # X AXIS
    # ========================================================

    ax.set_xlabel(
        "Rigidity [GV]",
        fontsize=30
    )

    ax.tick_params(
        axis="x",
        labelsize=24
    )

    # ========================================================
    # LEGEND
    # ========================================================

    allowed_patch = Patch(
        facecolor="white",
        edgecolor="none",
        label="Allowed"
    )

    forbidden_patch = Patch(
        facecolor="black",
        edgecolor="black",
        label="Forbidden"
    )

    ru_label = Line2D(
        [],
        [],
        color="none",
        label=f"Ru: {Ru:.3f} [GV]"
    )

    rc_label = Line2D(
        [],
        [],
        color="none",
        label=f"Rc: {Rc:.3f} [GV]"
    )

    rl_label = Line2D(
        [],
        [],
        color="none",
        label=f"Rl:  {Rl:.3f} [GV]"
    )

    ax.legend(
        handles=[
            allowed_patch,
            forbidden_patch,
            ru_label,
            rc_label,
            rl_label
        ],
        loc="upper right",
        fontsize=17,
        frameon=True,
        handletextpad=0.6
    )

    # ========================================================
    # BORDER
    # ========================================================

    for spine in ax.spines.values():

        spine.set_linewidth(
            1.5
        )

        spine.set_color(
            "black"
        )

    # ========================================================
    # LAYOUT
    # ========================================================

    plt.tight_layout()

    return fig, ax


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    stations_list = [
        "OULU"
    ]

    # ========================================================
    # RUN OTSO
    # ========================================================

    cone_results = cone(
        Stations=stations_list,
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params=rigidity_inputs,
        transmission_params=transmission_inputs,
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

    # ========================================================
    # EXTRACT RESULTS
    # ========================================================

    cone_df = cone_results[0]
    cutoff_df = cone_results[1]

    # ========================================================
    # PRINT RESULTS
    # ========================================================

    print()
    print("=" * 70)
    print("                         CONE RESULTS")
    print("=" * 70)

    print(
        cone_df
    )

    print()
    print("=" * 70)
    print("                     CUTOFF RESULTS")
    print("=" * 70)

    print(
        cutoff_df
    )

    # ========================================================
    # PLOT EACH STATION
    # ========================================================

    for station in stations_list:

        fig, ax = plot_cone(
            cone_df=cone_df,
            cutoff_df=cutoff_df,
            station=station,
            xmin=0.30,
            xmax=0.70
        )

        plt.show()