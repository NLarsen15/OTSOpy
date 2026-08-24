import json
import pandas as pd
import multiprocessing as mp
from collections import defaultdict

from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..data_classes.skymap_data import SkymapData
from ..utils import cutoff_utils as cu
from ..utils import cpus_utils as cpu_util


def FortranSkymap(Data: list, DataSkymap: SkymapData, queue: mp.Queue, cpus, JsonFile) -> None:

    cpu_util.set_process_affinity(cpus)

    # Group data by station name
    station_groups = defaultdict(list)

    for x in Data:
        station_groups[x[0]].append(x)

    if DataSkymap.model[1] == 99:
        mhd_utils.MHDinitialise(DataSkymap.MHDfile)

    # Dictionary of DataFrames
    all_results = {}

    for Station, station_data in station_groups.items():

        station_records = []

        for x in station_data:

            FortranData = fortran_data_utils.prepare_fortran_skymap(
                x,
                DataSkymap
            )

            returnvalues = cu.standard_cutoff(
                FortranData,
                DataSkymap.g,
                DataSkymap.h,
                DataSkymap
            )

            R_high, R_eff, R_low, transparency = returnvalues

            record = {
                "Zenith": float(x[4]),
                "Azimuth": float(x[5]),
                "Ru [GV]": float(R_high),
                "Rc [GV]": float(R_eff),
                "Rl [GV]": float(R_low),
                "PTF": float(transparency)
            }

            station_records.append(record)

            # update progress for every ZA calculation
        queue.put(1)

        all_results[Station] = pd.DataFrame(station_records)

    json_results = {
        station: df.to_dict(orient="records")
        for station, df in all_results.items()
    }

    with open(JsonFile, "w", encoding="utf-8") as f:
        json.dump(
            json_results,
            f,
            indent=4
        )

    return