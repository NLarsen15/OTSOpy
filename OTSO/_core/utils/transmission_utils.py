import pandas as pd
import numpy as np
import os

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils import cutoff_utils as cu

def transmission(FortranData, R_high, R_low, Station, Data):

    full_rigidities = np.arange(
        Data.endrigidity + Data.rigiditystep,
        Data.startrigidity + Data.rigiditystep,
        Data.rigiditystep,
    )

    decimals = len(str(FortranData.rigiditystep).split(".")[1])

    prefilled_df = pd.DataFrame({"R [GV]": full_rigidities})

    prefilled_df[f"{Station}_TF"] = np.where(
        prefilled_df["R [GV]"] < R_low,
        0.0,
        np.where(prefilled_df["R [GV]"] > R_high, 1.0, np.nan)
    )

    FortranData.transmissionrres = Data.transmissionRstep
    FortranData.transmissionsamples = Data.transmissionsamples

    FortranData.startrigidity = R_high + FortranData.rigiditystep
    FortranData.endrigidity = R_low - FortranData.rigiditystep

    if FortranData.endrigidity < 0.0:
        FortranData.endrigidity = 0.0

    FortranData.n = cu.compute_number_of_rigidities(FortranData.startrigidity, 
                                                    FortranData.endrigidity, FortranData.rigiditystep)

    g = Data.g
    h = Data.h

    Transmissions = np.zeros(FortranData.n, dtype=np.float64)
    Rigidities = np.zeros(FortranData.n, dtype=np.float64)

    OTSOLib.transmission(
        FortranData,
        g,
        h,
        Rigidities,
        Transmissions
    )

    transmission_df = pd.DataFrame({
      "R [GV]": Rigidities,
      f"{Station}_TF": Transmissions})

    transmission_df = transmission_df.sort_values(by="R [GV]").reset_index(drop=True)

    calculated_df = pd.DataFrame({
        "R [GV]": Rigidities,
        f"{Station}_TF": Transmissions
    })

    prefilled_df["R [GV]"] = prefilled_df["R [GV]"].round(decimals)
    calculated_df["R [GV]"] = calculated_df["R [GV]"].round(decimals)

    Transmissiondf = prefilled_df.set_index("R [GV]")

    Transmissiondf.update(calculated_df.set_index("R [GV]"))

    Transmissiondf = (
        Transmissiondf
        .reset_index()
        .sort_values("R [GV]")
        .reset_index(drop=True)
    )

    Transmissiondf.update(calculated_df.set_index("R [GV]"))

    Transmissiondf = Transmissiondf.sort_values("R [GV]")

    return Transmissiondf

def transmission_full(FortranData, R_high, R_low, Station, Data):

    FortranData.transmissionrres = Data.transmissionRstep
    FortranData.transmissionsamples = Data.transmissionsamples

    if FortranData.endrigidity < 0.0:
        FortranData.endrigidity = 0.0

    FortranData.n = cu.compute_number_of_rigidities(
        FortranData.startrigidity,
        FortranData.endrigidity,
        FortranData.rigiditystep
    )

    g = Data.g
    h = Data.h

    Transmissions = np.zeros(FortranData.n, dtype=np.float64)
    Rigidities = np.zeros(FortranData.n, dtype=np.float64)

    OTSOLib.transmission(
        FortranData,
        g,
        h,
        Rigidities,
        Transmissions
    )

    decimals = len(str(FortranData.rigiditystep).split(".")[1])

    Transmissiondf = pd.DataFrame({
        "R [GV]": np.round(Rigidities, decimals),
        f"{Station}_TF": Transmissions
    })

    Transmissiondf = (
        Transmissiondf
        .sort_values("R [GV]")
        .reset_index(drop=True)
    )

    return Transmissiondf
