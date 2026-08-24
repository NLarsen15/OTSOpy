import pandas as pd
import numpy as np
import os

from ..libs.MiddleMan import Middleman as OTSOLib


def find_cutoff_range(df, rigidity_step):
    # Ensure sorted by rigidity
    df = df.sort_values("R [GV]").reset_index(drop=True)
    R_step = rigidity_step
    decimals = len(str(R_step).split(".")[1])


    # Lowest rigidity where A == 1
    R_low = df.loc[df["Filter"] == 1, "R [GV]"].min()

    # Highest forbidden rigidity (A < 0)
    R_last_forbidden = df.loc[df["Filter"] < 0, "R [GV]"].max()

    # First rigidity above the last forbidden value
    R_high = df.loc[df["R [GV]"] > R_last_forbidden, "R [GV]"].min()

    # Penumbra region (exclude R_low and R_high)
    penumbra = df[(df["R [GV]"] >= R_low) & (df["R [GV]"] < R_high)]

    # Number of allowed values in the penumbra
    count_allowed = (penumbra["Filter"] == 1).sum()

    # Total number of rigidity values in the penumbra
    count_total = len(penumbra)

    # Penumbral transparency
    transparency = count_allowed / count_total if count_total > 0 else float("nan")
    if np.isnan(transparency):
        transparency = 0.0

    # Effective cutoff
    R_eff = R_high - R_step * count_allowed

    R_eff = round(R_eff, decimals)
    R_low = round(R_low, decimals)
    R_high = round(R_high, decimals)
    transparency = round(transparency, decimals)


    return (
        R_low,
        R_high,
        transparency,
        R_eff
    )

def compute_number_of_rigidities(start_rigidity, end_rigidity, rigidity_step):

    #print(f"Computing number of rigidities with start: {start_rigidity}, end: {end_rigidity}, step: {rigidity_step}")
    
    n = int((start_rigidity - end_rigidity) / rigidity_step)

    return n

def Rigidity_scan(FortranData, g, h, CutoffDataInstance):
    smallscan = False
    FortranData.rigiditystep = 0.25

    #print("scan")
    #print(
    #    "PID:", os.getpid(),
    #    "CPUs:", os.sched_getaffinity(0)
    #)

    n = compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n
    FortranData.gyropercent = CutoffDataInstance.maxsteppercent

    Rigidities = np.zeros(n, dtype=np.float64)
    Allowed = np.zeros(n, dtype=np.int32)
    
    # Call new Fortran interface

    OTSOLib.cutoff(
        FortranData,
        g,
        h,
        Rigidities,
        Allowed
    )

    df = pd.DataFrame({
    "R [GV]": Rigidities,
    "Filter": Allowed})

    if np.all(Allowed == 1):
        smallscan = smallest_step_check(FortranData, g, h, CutoffDataInstance)
        if smallscan:
           return FortranData, smallscan

  # Sort by rigidity
    df = df.sort_values(by="R [GV]").reset_index(drop=True)

    R_low, R_high, transparency, R_eff = find_cutoff_range(df, FortranData.rigiditystep)

    
    # Create result DataFrame with rigidities

    if np.isnan(R_low):
        R_low = 0
    if np.isnan(R_high):
        R_high = R_low

    R_low = R_low - 0.75
    if (R_low < 0):
       R_low = 0

    R_high = R_high + 0.75
    if (R_high > FortranData.startrigidity):
       R_high = FortranData.startrigidity

    FortranData.rigiditystep = CutoffDataInstance.rigiditystep
    FortranData.startrigidity = R_high
    FortranData.endrigidity = R_low

    n = compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n

    return FortranData, smallscan

def standard_cutoff(FortranData, g, h, Data):
    smallscan = False

    if Data.rigidityscan == "ON":
        FortranData, smallscan = Rigidity_scan(FortranData, g, h, Data)

    if smallscan:
            returnvalues = [0.0, 0.0, 0.0, 0.0]
            R_low, R_high, transparency, R_eff = returnvalues
            
    else:
            Rigidities = np.zeros(FortranData.n, dtype=np.float64)
            Allowed = np.zeros(FortranData.n, dtype=np.int32)
            FortranData.gyropercent = Data.maxsteppercent
            
            # Call new Fortran interface
            OTSOLib.cutoff(
                FortranData,
                g,
                h,
                Rigidities,
                Allowed
            )

            df = pd.DataFrame({
            "R [GV]": Rigidities,
            "Filter": Allowed})

            # Sort by rigidity
            df = df.sort_values(by="R [GV]").reset_index(drop=True)

            R_low, R_high, transparency, R_eff = find_cutoff_range(df, FortranData.rigiditystep)

    returnvalues = [R_high, R_eff, R_low, transparency]

    return returnvalues

def Apparent_cutoff(Positioninfo, FortranData, g, h, Data):

    Zenith_and_Azimuths = [
        [0, 0],        # vertical
        [30, 0],
        [30, 45],
        [30, 90],
        [30, 135],
        [30, 180],
        [30, 225],
        [30, 270],
        [30, 315]
    ]

    # Store results from each direction
    results = []

    for za in Zenith_and_Azimuths:
        smallscan = False
        Position = [Positioninfo[0], Positioninfo[1], Positioninfo[2], za[0], za[1]]
        FortranData.positionin = Position
        FortranData.startrigidity = Data.startrigidity
        FortranData.endrigidity = Data.endrigidity

        if Data.rigidityscan == "ON":
            FortranData, smallscan = Rigidity_scan(FortranData, Data)

        if smallscan:
            R_high, R_eff, R_low, transparency = [0.0, 0.0, 0.0, 0.0]
        else:
            n = compute_number_of_rigidities(
                FortranData.startrigidity,
                FortranData.endrigidity,
                FortranData.rigiditystep
            )

            FortranData.n = n
            FortranData.gyropercent = Data.maxsteppercent

            Rigidities = np.zeros(n, dtype=np.float64)
            Allowed = np.zeros(n, dtype=np.int32)

            OTSOLib.cutoff(
                FortranData,
                g,
                h,
                Rigidities,
                Allowed
            )

            df = pd.DataFrame({
                "R [GV]": Rigidities,
                "Filter": Allowed
            }).sort_values(by="R [GV]").reset_index(drop=True)

            R_low, R_high, transparency, R_eff = find_cutoff_range(df, FortranData.rigiditystep)

        results.append({
            "Direction": tuple(za),
            "Ru": R_high,
            "Rc": R_eff,
            "Rl": R_low,
            "TF": transparency
            })

    rigidity_dataframe = pd.DataFrame(results).set_index("Direction")

    # Weights: vertical = 1/2, each of the other 8 = 1/16
    weights = np.array([0.5] + [1/16] * 8)

    R_step = FortranData.rigiditystep
    decimals = len(str(R_step).split(".")[1])

    weighted_mean = (
        rigidity_dataframe[["Ru", "Rc", "Rl", "TF"]]
        .multiply(weights, axis=0)
        .sum(axis=0)
        .round(decimals)
    )

    returnvalues = [weighted_mean["Ru"], weighted_mean["Rc"], 
                    weighted_mean["Rl"], weighted_mean["TF"]]

    return returnvalues


def smallest_step_check(FortranData, g, h, CutoffDataInstance):
    n = 1

    FortranData.gyropercent = CutoffDataInstance.maxsteppercent

    FortranData.n = n
    FortranData.rigiditystep = CutoffDataInstance.rigidityarray[2] + CutoffDataInstance.rigidityarray[1]
    #FortranData.rigiditystep = 0.05

    Rigidities = np.zeros(n, dtype=np.float64)
    Allowed = np.zeros(n, dtype=np.int32)
    
    # Call new Fortran interface
    OTSOLib.cutoff(
        FortranData,
        g,
        h,
        Rigidities,
        Allowed
    )

    df = pd.DataFrame({
    "R [GV]": Rigidities,
    "Filter": Allowed})

    if Allowed[0] == 1:
       result = True # The lowest tested rigidity is allowed
       #print("Smallest Step Escaped")
    else:
       result = False # The lowest tested rigidity is not allowed
       #print("Smallest Step Did Not Escape")

    return result