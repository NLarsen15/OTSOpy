import pandas as pd
import multiprocessing as mp
import numpy as np
import json

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..data_classes.cone_data import ConeData
from ..utils import cutoff_utils as cu
from ..utils import transmission_utils as tu
from ..utils import cpus_utils as cpu_util


def FortranCone(Data: list, ConeDataInstance: ConeData, queue: mp.Queue, cpus, JsonFile, lock) -> None:

    cpu_util.set_process_affinity(cpus)

    if ConeDataInstance.model[1] == 99:
        mhd_utils.MHDinitialise(ConeDataInstance.MHDfile)

    for x in Data:

        Station = x[0]

        FortranData = fortran_data_utils.prepare_fortran_cone(
            x, ConeDataInstance
        )

        g = ConeDataInstance.g
        h = ConeDataInstance.h

        Rigidities = np.zeros(FortranData.n, dtype=np.float64)
        Allowed = np.zeros(FortranData.n, dtype=np.int32)
        Asymlat = np.zeros(FortranData.n, dtype=np.float64)
        Asymlong = np.zeros(FortranData.n, dtype=np.float64)

        OTSOLib.cone(
            FortranData,
            g,
            h,
            Rigidities,
            Allowed,
            Asymlat,
            Asymlong,
        )

        decimals = len(str(FortranData.rigiditystep).split(".")[1])

        Conedf = (
            pd.DataFrame(
                {
                    "R [GV]": Rigidities.round(decimals),
                    "Filter": Allowed,
                    "Alat": Asymlat,
                    "Along": Asymlong,
                }
            )
            .sort_values("R [GV]")
            .reset_index(drop=True)
        )

        Conedf["Alat"] = Conedf["Alat"].round(4)
        Conedf["Along"] = Conedf["Along"].round(4)

        R_low, R_high, transparency, R_eff = cu.find_cutoff_range(
            Conedf, FortranData.rigiditystep
        )

        returnvalues = [R_high, R_eff, R_low, transparency]

        Rigiditydataframe = pd.DataFrame(
            {Station: returnvalues},
            index=["Ru", "Rc", "Rl", "PTF"],
        )

        

        # ----------------------------------------------------------
        # Append transmission (if requested)
        # ----------------------------------------------------------
        if ConeDataInstance.transmission:

            Transmissiondf = tu.transmission(
                FortranData,
                R_high,
                R_low,
                Station,
                ConeDataInstance,
            )

            Conedf = Conedf.merge(
                Transmissiondf,
                on="R [GV]",
                how="left",
            )

            tf_col = next(
                col for col in Transmissiondf.columns if col != "R [GV]"
            )

            Conedf[Station] = (
                Conedf["Filter"].astype(str)
                + ConeDataInstance.delim
                + Conedf["Alat"].map("{:.4f}".format)
                + ConeDataInstance.delim
                + Conedf["Along"].map("{:.4f}".format)
                + ConeDataInstance.delim
                + Conedf[tf_col].astype(str)
            )

        else:

            Conedf[Station] = (
                  Conedf["Filter"].astype(str)
                  + ConeDataInstance.delim
                  + Conedf["Alat"].map("{:.4f}".format)
                  + ConeDataInstance.delim
                  + Conedf["Along"].map("{:.4f}".format)
              )

        Conedf = Conedf[["R [GV]", Station]]

        record = { 
              "station": Station,
              "cone": Conedf.to_dict(orient="split"),
              "rigidities": Rigiditydataframe.to_dict(orient="split"),
              }

        with lock:
          with open(JsonFile, "a", encoding="utf-8") as f:
              json.dump(record, f)
              f.write("\n")

        queue.put(1)

    return


