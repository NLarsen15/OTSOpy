from OTSO import transmission

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import StepPatch


# ============================================================
# 1. INPUT PARAMETERS
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
    "rigiditystep": 0.01
}

transmission_inputs = {
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
    "corenum": 5,
    "threadnum": 1,
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
    "inputcoord": "GDZ"
}

data_retrieval_inputs = {
    "serverdata": "OFF",
    "livedata": "OFF"
}


# ============================================================
# 2. STATIONS / CUSTOM LOCATIONS
# ============================================================

stations_list = [
    "ROME"
]


# ============================================================
# 3. RUN OTSO
# ============================================================

print("Running OTSO transmission calculation...")

transmission_results = transmission(
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

print("OTSO calculation finished.")


# ============================================================
# 4. GET TRANSMISSION DATAFRAME
# ============================================================

# The first result is the transmission DataFrame
df = transmission_results[0]

print("\nColumns returned by OTSO:")
print(df.columns.tolist())


# ============================================================
# 5. FIND THE TRANSMISSION FUNCTION COLUMN
# ============================================================

# OTSO names the transmission column according to the station,
# e.g. ROME_TF, OULU_TF, ATHN_TF, CALG_TF.

tf_columns = [
    column for column in df.columns
    if column.endswith("_TF")
]

if len(tf_columns) == 0:
    raise ValueError(
        "No transmission-function column (*_TF) was found.\n"
        f"Columns returned by OTSO: {df.columns.tolist()}"
    )

# If only one station is being used, this is the one we want.
tf_column = tf_columns[0]

print(f"\nUsing transmission column: {tf_column}")


# ============================================================
# 6. EXTRACT RIGIDITY AND TRANSMISSION FUNCTION
# ============================================================

R = df["R [GV]"].to_numpy()
TF = df[tf_column].to_numpy()


# Remove NaN values if present
valid = np.isfinite(R) & np.isfinite(TF)

R = R[valid]
TF = TF[valid]


# Sort by rigidity
sort_index = np.argsort(R)

R = R[sort_index]
TF = TF[sort_index]


# ============================================================
# 7. CALCULATE STEP-PLOT BIN EDGES
# ============================================================

dR = np.median(np.diff(R))

edges = np.concatenate(
    (
        [R[0] - dR / 2],
        R + dR / 2
    )
)


# ============================================================
# 8. CREATE PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(12, 5))

patch = StepPatch(
    values=TF,
    edges=edges,
    baseline=0,
    facecolor="white",
    edgecolor="black",
    hatch="///",
    linewidth=1.5
)

ax.add_patch(patch)


# ============================================================
# 9. PLOT FORMATTING
# ============================================================

ax.set_xlim(5.5, 6.75)
ax.set_ylim(0, 1.0)

ax.set_xlabel(
    "Rigidity [GV]",
    fontsize=18
)

ax.set_ylabel(
    "Transmission Function",
    fontsize=18
)

ax.tick_params(
    axis="both",
    labelsize=16
)

plt.tight_layout()


# ============================================================
# 10. SAVE FIGURE
# ============================================================

plt.savefig(
    "transmission.png",
    dpi=400,
    bbox_inches="tight"
)

plt.show()