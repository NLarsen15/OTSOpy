import json
import pandas as pd
import multiprocessing as mp
import numpy as np

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..data_classes.planet_data import PlanetData
from ..utils import cutoff_utils as cu
from ..utils import transmission_utils as tu
from ..utils import cpus_utils as cpu_util


def FortranPlanet(
    TaskQueue: mp.Queue,
    DataPlanet: PlanetData,
    queue: mp.Queue,
    cpus,
    JsonFile
) -> None:

    cpu_util.set_process_affinity(cpus)

    if DataPlanet.model[1] == 99:
        mhd_utils.MHDinitialise(DataPlanet.MHDfile)

    with open(JsonFile, "a", encoding="utf-8") as f:

        while True:

            x = TaskQueue.get()

            if x is None:
                break

            Station = f"{x[1]}_{x[2]}"
            Position = [
                x[3],
                x[1],
                x[2],
                x[4],
                x[5]
            ]

            # --------------------------------------------------
            # Prepare Fortran data
            # --------------------------------------------------

            FortranData = (
                fortran_data_utils.prepare_fortran_cutoff(
                    x,
                    DataPlanet
                )
            )

            # --------------------------------------------------
            # Cutoff calculation
            # --------------------------------------------------

            if DataPlanet.cutoff_comp == "Apparent":

                returnvalues = cu.Apparent_cutoff(
                    Position,
                    FortranData,
                    DataPlanet.g,
                    DataPlanet.h,
                    DataPlanet
                )

            else:

                returnvalues = cu.standard_cutoff(
                    FortranData,
                    DataPlanet.g,
                    DataPlanet.h,
                    DataPlanet
                )

            R_high, R_eff, R_low, transparency = returnvalues

            # --------------------------------------------------
            # Transmission
            # --------------------------------------------------

            Transmissiondf = None

            if DataPlanet.transmission:

                Transmissiondf = tu.transmission(
                    FortranData,
                    R_high,
                    R_low,
                    Station,
                    DataPlanet
                )

            # --------------------------------------------------
            # Asymptotic calculation
            # --------------------------------------------------

            AsyAllowed = None
            AsyLat = None
            AsyLong = None

            if DataPlanet.asymptotic == "YES":

                Energy_List = DataPlanet.asymlevels.copy()
                P_List = []

                if DataPlanet.unit == "GeV":

                    E_0 = 0.938

                    for i in Energy_List:

                        R = (i**2 + 2 * i * E_0) ** 0.5
                        P_List.append(R)

                elif DataPlanet.unit == "GV":

                    P_List = Energy_List.copy()

                AsyAllowed = np.zeros(
                    len(P_List),
                    dtype=np.int32
                )

                AsyLat = np.zeros(
                    len(P_List),
                    dtype=np.float64
                )

                AsyLong = np.zeros(
                    len(P_List),
                    dtype=np.float64
                )

                Plist_array = np.array(
                    P_List,
                    dtype=np.float64
                )

                OTSOLib.trajectory(
                    FortranData,
                    DataPlanet.g,
                    DataPlanet.h,
                    Plist_array,
                    len(P_List),
                    AsyAllowed,
                    AsyLat,
                    AsyLong
                )

            # --------------------------------------------------
            # Create result record
            # --------------------------------------------------

            record = {
                "Latitude": x[1],
                "Longitude": x[2],
                "Ru [GV]": R_high,
                "Rc [GV]": R_eff,
                "Rl [GV]": R_low,
                "PTF": transparency
            }

            # --------------------------------------------------
            # Add asymptotic results
            # --------------------------------------------------

            if DataPlanet.asymptotic == "YES":

                asymlevels_with_units = [
                    f"{level} [{DataPlanet.unit}]"
                    for level in DataPlanet.asymlevels
                ]

                formatted_asymptotic_data = [
                    f"{AsyAllowed[i]}{DataPlanet.delim}"
                    f"{round(AsyLat[i], 4)}{DataPlanet.delim}"
                    f"{round(AsyLong[i], 4)}"
                    for i in range(len(AsyAllowed))
                ]

                asymptotic_df = pd.DataFrame(
                    [[Station] + formatted_asymptotic_data],
                    columns=[
                        "Station"
                    ] + asymlevels_with_units
                )

                record["asymptotic"] = (
                    asymptotic_df.to_dict(
                        orient="records"
                    )
                )

            # --------------------------------------------------
            # Add transmission results
            # --------------------------------------------------

            if DataPlanet.transmission:

                transmission = {}

                tf_column = Transmissiondf.columns[1]

                for _, row in Transmissiondf.iterrows():

                    transmission[
                        f'{row["R [GV]"]} [GV]'
                    ] = row[tf_column]

                record["transmission"] = transmission

            # --------------------------------------------------
            # Write JSON record
            # --------------------------------------------------

            json.dump(record, f)
            f.write("\n")

            queue.put(1)


